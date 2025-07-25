restoredefaultpath;

clear all;
clc;
close all;
warning off;

%primary assumptions 

hours_shift = [0, 4, -4];

mu = 0.85; % roundtrip efficiency
power_to_energy_ratio = 0.25; 

% battery limits
SOE_min = 0.1;
SOE_max = 0.9;
starting_SoE = 0.5;

battery_present = true;

pack_sizes = [0 100.0];

% load in data
capacity_profiles_raw = readtable("../data/gac_sol_1.0mw_capacity_profiles.csv");
data_table = readtable("../data/excess_loadDemand_emissions_data.xlsx", "Sheet", 'Sheet1');

% process data
len_hours_shift = length(hours_shift);
len_pack_sizes = length(pack_sizes);

pts_capacity = circshift(capacity_profiles_raw{:,6:end}, -5, 2);
[table_len, N] = size(pts_capacity);

days_of_year = (1:365)';
sat_days = days_of_year(6:7:end);
sun_days = days_of_year(7:7:end);

sat_hrs_idxs = return_hours_of_day(sat_days);
sun_hrs_idxs = return_hours_of_day(sun_days);

sat_hrs_idxs = sat_hrs_idxs(:);
sun_hrs_idxs = sun_hrs_idxs(:);

ell_scaling_factor = 1.0; % scaling factor for aggressive electrification scenario
% ell_scaling_factor = 801/(848 + 460 + 801); % scaling factor for VMT-based electrification scenario
ell_t_base = (ell_scaling_factor.*data_table.Grid_Energy_hrs_kWh_(1:N))./1e3; % grid demand in megawatts
ell_t_base(sat_hrs_idxs) = 0.0;
ell_t_base(sun_hrs_idxs) = 0.0;

co2_t_base = (ell_scaling_factor.*data_table.Grid_CO2_hrs_kg_(1:N))./907.185; % CO2 grid production in US ton
co2_to_ell_base(ell_t_base == 0.0) = 0.0;

% set up
site_pts_id = (1:table_len)';
hours_in_year = (1:N)';

cycles_optimal = zeros(len_pack_sizes, table_len);

tot_solar = zeros(table_len, 1);
total_load = zeros(table_len, 1);
co2_emissions_baseline = zeros(table_len, 1);

carbon_subtracted_optimal = zeros(len_pack_sizes, table_len);
carbon_reduction_optimal = zeros(len_pack_sizes, table_len);

tot_solar_curtailed_optimal = zeros(len_pack_sizes, table_len);
tot_solar_to_batt_optimal = zeros(len_pack_sizes, table_len);
tot_batt_dchg_optimal = zeros(len_pack_sizes, table_len);
tot_batt_chg_optimal = zeros(len_pack_sizes, table_len);
tot_renewables_optimal = zeros(len_pack_sizes, table_len);

solar_fractions_optimal = zeros(len_pack_sizes, table_len);
ESS_fractions_optimal = zeros(len_pack_sizes, table_len);
grid_fractions_optimal = zeros(len_pack_sizes, table_len);

solar_fractions_ESSCharge_optimal = zeros(len_pack_sizes, table_len);
grid_fractions_ESSCharge_optimal = zeros(len_pack_sizes, table_len);

zeros_vec = zeros(N, 1); % for convenience later

% main loop
variable_names = cell(13, 1);
variable_names(:) = {'CO2_intensity_t'; ...
    'Pdemand_total_MW'; ...
    'Ppv_total_MW'; ...
    'Pgrid_total_MW'; ...
    'Pgrid_to_load_MW'; ...
    'Pgrid_to_ESS_MW'; ...
    'Ppv_to_load_MW'; ...
    'Ppv_to_ESS_MW'; ...
    'Ppv_crtail_MW'; ...
    'Pess_chg_MW'; ...
    'Pess_dchg_MW'; ...
    'Eess_MWh'; ...
    'dchg_indicator'};

summary_variable_names = cell(16, 1);
summary_variable_names(:) = {'sc_gid'; ...
    'total_load_MWh'; ...
    'total_CO2_baseline_t'; ...
    'CO2_removed_pcent'; ...
    'CO2_removed_t'; ...
    'solar_fraction_num'; ...
    'ESS_fraction_num'; ...
    'grid_fraction_num'; ...
    'solar_chargeFraction_num'; ...
    'grid_chargeFraction_num'; ...
    'tot_solar_MWh'; ...
    'tot_solar_to_batt_MWh'; ...
    'tot_solar_curtailed_MWh'; ...
    'tot_batt_dchg_MWh'; ...
    'tot_batt_chg_MWh'; ...
    'ESS_cycles_num'};

opts = optimoptions('intlinprog', ...
                    'Display', 'none', ...
                    'Heuristics', 'advanced');

len_pack_sizes = length(pack_sizes);

% create initial solutions - assume no microgrid; use only grid
x0(len_pack_sizes, table_len, len_hours_shift) = struct();
for s=1:len_hours_shift

    hours_to_shift = hours_shift(s);

    ell_t = circshift(ell_t_base(1:N), hours_to_shift);

    for k=1:table_len

        pv_cap_t_base = pts_capacity(k, :);
        pv_t = reshape(pv_cap_t_base(1:N), [], 1);

        for ps=1:len_pack_sizes
            pack_size = pack_sizes(ps); % size of pack in MWh

            % % trial initial solution - no battery
            x0(ps, k, s).g_t = ell_t;
            x0(ps, k, s).g_load = ell_t;
            x0(ps, k, s).pv_crtail = pv_t;
            x0(ps, k, s).E_t = starting_SoE.*pack_size.*battery_present.*ones(N+1, 1);

            x0(ps, k, s).g_ESS = zeros_vec;
            x0(ps, k, s).pv_load = zeros_vec;
            x0(ps, k, s).pv_ESS = zeros_vec;
            x0(ps, k, s).c_t = zeros_vec;
            x0(ps, k, s).d_t = zeros_vec;
            x0(ps, k, s).dbin_t = zeros(N, 1);
        end
    end
end

for s=1:len_hours_shift

    hours_to_shift = hours_shift(s);

    ell_t = circshift(ell_t_base(1:N), hours_to_shift);
    co2_t = circshift(co2_t_base(1:N), hours_to_shift);
    
    nonzero_ell_t_idxs = find(ell_t ~= 0.0);
    zero_ell_t_idxs = (ell_t == 0.0);
    
    out_fname_pre = ['../data/results/site_varied/full_output_shift_' num2str(hours_to_shift, '%.0f') '_mu_' num2str(mu, "%.2f") '_P2Eratio_' num2str(power_to_energy_ratio, '%.2f')];
    summary_fpath = ['../data/results/site_varied/summary_shift_' num2str(hours_to_shift, '%.0f') '_mu_' num2str(mu, "%.2f") '_P2Eratio_' num2str(power_to_energy_ratio, '%.2f') '_outmat.xlsx'];
    
    parfor k=1:table_len
    
        pv_cap_t_base = pts_capacity(k, :);
        pv_t = reshape(pv_cap_t_base(1:N), [], 1);
    
        tot_solar(k) = sum(pv_t);
    
        total_load(k) = sum(ell_t);
        co2_emissions_baseline(k) = sum(co2_t);
    
        filename = [out_fname_pre '_SolarFarmSite_' num2str(k, '%03.0f') '_outmat.xlsx']; % output file name
    
        for ps=1:len_pack_sizes
    
            pack_size = pack_sizes(ps);

            P_ESS = power_to_energy_ratio.*pack_size; % maximum power in MW
    
            % ============= OPTIMAL PROBLEM =============
    
            % initialize the problem
            prob = optimproblem('ObjectiveSense', 'maximize');
    
            % decision variables
            g_t       = optimvar('g_t',       N,   1, 'Type', 'continuous', 'LowerBound', 0.0);
            g_load    = optimvar('g_load',    N,   1, 'Type', 'continuous', 'LowerBound', 0.0);
            g_ESS     = optimvar('g_ESS',     N,   1, 'Type', 'continuous', 'LowerBound', 0.0);
            pv_load   = optimvar('pv_load',   N,   1, 'Type', 'continuous', 'LowerBound', 0.0);
            pv_ESS    = optimvar('pv_ESS',    N,   1, 'Type', 'continuous', 'LowerBound', 0.0);
            pv_crtail = optimvar('pv_crtail', N,   1, 'Type', 'continuous', 'LowerBound', 0.0);
            c_t       = optimvar('c_t',       N,   1, 'Type', 'continuous', 'LowerBound', 0.0);
            d_t       = optimvar('d_t',       N,   1, 'Type', 'continuous', 'LowerBound', 0.0);
            E_t       = optimvar('E_t',       N+1, 1, 'Type', 'continuous', 'LowerBound', SOE_min.*pack_size.*battery_present, 'UpperBound', SOE_max.*pack_size.*battery_present);
            dbin_t    = optimvar('dbin_t',    N,   1, 'Type', 'integer',    'LowerBound', 0,                                   'UpperBound', 1);
    
            define objective function
            prob.Objective = 1e2.*(1.0 - (...
                sum(g_t(nonzero_ell_t_idxs).*co2_t(nonzero_ell_t_idxs)./ell_t(nonzero_ell_t_idxs)) ...
                + 0.048.*sum(pv_load + pv_ESS) ...
                )./co2_emissions_baseline(k));

            % inequality constraints
            prob.Constraints.ineq1 = d_t <= dbin_t.*P_ESS.*battery_present;
            prob.Constraints.ineq2 = c_t <= (1 - dbin_t).*P_ESS.*battery_present;
    
            % equality constraints
            prob.Constraints.eq1 = ell_t                == g_load + pv_load + d_t;
            prob.Constraints.eq2 = pv_t                 == pv_load + pv_ESS + pv_crtail;
            prob.Constraints.eq3 = g_t                  == g_ESS + g_load;
            prob.Constraints.eq4 = c_t                  == pv_ESS + g_ESS;
            prob.Constraints.eq5 = E_t(2:end)           == E_t(1:end-1) + sqrt(mu).*c_t - d_t./sqrt(mu);
            prob.Constraints.eq6 = E_t([1; end])        == starting_SoE.*pack_size.*battery_present;
            prob.Constraints.eq7 = g_t(zero_ell_t_idxs) == 0.0;
    
            % solve the problem
            [sol,fval] = solve(prob, x0(ps, k, s), 'Solver', 'intlinprog', 'Options', opts);
    
            % write solution to excel file
            output_array = [co2_t ell_t pv_t sol.g_t sol.g_load sol.g_ESS sol.pv_load sol.pv_ESS sol.pv_crtail sol.c_t sol.d_t sol.E_t(2:end) sol.dbin_t];
            output_table = array2table(output_array, 'VariableNames', variable_names);
            writetable(output_table, filename, 'Sheet', [num2str(pack_size, '%.0f') ' MWh Battery - Optimal'], 'WriteVariableNames', true);
    
            % Post-processing
            cycles_optimal(ps, k) = compute_cycles([0.0; hours_in_year], sol.E_t./pack_size);

            carbon_subtracted_optimal(ps, k) = co2_emissions_baseline(k) - (sum(sol.g_t(nonzero_ell_t_idxs).*co2_t(nonzero_ell_t_idxs)./ell_t(nonzero_ell_t_idxs)) + 0.048.*sum(sol.pv_load + sol.pv_ESS));
            carbon_reduction_optimal(ps, k) = fval;
    
            tot_solar_curtailed_optimal(ps, k) = sum(sol.pv_crtail);
            tot_solar_to_batt_optimal(ps, k) = sum(sol.c_t);
            tot_batt_dchg_optimal(ps, k) = sum(sol.d_t);
            tot_batt_chg_optimal(ps, k) = sum(sol.c_t);
    
            solar_fractions_optimal(ps, k) = sum(sol.pv_load(nonzero_ell_t_idxs)./ell_t(nonzero_ell_t_idxs))./length(nonzero_ell_t_idxs);
            ESS_fractions_optimal(ps, k) = sum(sol.d_t(nonzero_ell_t_idxs)./ell_t(nonzero_ell_t_idxs))./length(nonzero_ell_t_idxs);
            grid_fractions_optimal(ps, k) = sum(sol.g_load(nonzero_ell_t_idxs)./ell_t(nonzero_ell_t_idxs))./length(nonzero_ell_t_idxs);
    
            nonzero_charge_idxs = find(sol.c_t ~= 0.0);
            solar_fractions_ESSCharge_optimal(ps, k) = sum(sol.pv_ESS(nonzero_charge_idxs)./sol.c_t(nonzero_charge_idxs))./length(nonzero_charge_idxs);
            grid_fractions_ESSCharge_optimal(ps, k) = sum(sol.g_ESS(nonzero_charge_idxs)./sol.c_t(nonzero_charge_idxs))./length(nonzero_charge_idxs);

        end
    end

    % write summary file
    for ps=1:length(pack_sizes)

        pack_size = pack_sizes(ps);
    
        optimal_summary_table = array2table([site_pts_id(:), ...
            total_load(:), ...
            co2_emissions_baseline(:), ...
            carbon_reduction_optimal(ps,:)', ...
            carbon_subtracted_optimal(ps,:)', ...
            solar_fractions_optimal(ps,:)', ...
            ESS_fractions_optimal(ps,:)', ...
            grid_fractions_optimal(ps,:)', ...
            solar_fractions_ESSCharge_optimal(ps,:)', ...
            grid_fractions_ESSCharge_optimal(ps,:)', ...
            tot_solar(:), ...
            tot_solar_to_batt_optimal(ps,:)', ...
            tot_solar_curtailed_optimal(ps,:)', ...
            tot_batt_dchg_optimal(ps,:)', ...
            tot_batt_chg_optimal(ps,:)', ...
            cycles_optimal(ps,:)'], 'VariableNames', summary_variable_names);
    
        writetable(optimal_summary_table, summary_fpath, 'Sheet', [num2str(pack_size, '%.0f') ' MWh Battery - Optimal'], 'WriteVariableNames', true);

    end

end

function hrs = return_hours_of_day(day_idxs)

    hrs = zeros(24, length(day_idxs));

    for i=1:length(day_idxs)

        n = day_idxs(i);

        hrs_low = 24.*(n-1) + 1;
        hrs_high = 24.*n;
        
        hrs(:,i) = (hrs_low:hrs_high)';
    end

end
