clear;
clc;
close all;
warning off;

%% user-set parameters

N = 8760; % number of hours in to consider in dispatch (max: 8760 h = 1yr)

hours_shift = [-4; 0; 4]; % demand shifts to consider

% battery system parameters
battery_present = true;

mu = 0.85; % roundtrip efficiency
power_to_energy_ratio = 0.25; 

SOE_min = 0.10;
SOE_max = 0.90;
starting_SoE = 0.5;

pack_sizes = [0.0; 100.0]; % battery system sizes (MW) to consider 

%% prepare arrays to be used later
hours_in_year = (1:N)';
days_of_year = (1:365)';
sat_days = days_of_year(6:7:end);
sun_days = days_of_year(7:7:end);

sat_hrs_idxs = return_hours_of_day(sat_days);
sun_hrs_idxs = return_hours_of_day(sun_days);

sat_hrs_idxs = sat_hrs_idxs(:);
sun_hrs_idxs = sun_hrs_idxs(:);

zeros_vec = zeros(N, 1); % for convenience

%% load necessary data

pv_data_table = readtable("gac_sol_1.0mw_capacity_profiles.csv");
[table_len, table_width] = size(pv_data_table);

data_table = readtable("Sample data for battery use optimization r3.xlsx", "Sheet", 'Median temp');

ell_t_base = 2.0.*data_table.Grid_Energy_hrs_kWh_(1:N)./1e3; % grid demand in megawatts
% ensure weekends are zeroed-out
ell_t_base(sat_hrs_idxs) = 0.0;
ell_t_base(sun_hrs_idxs) = 0.0;

co2_t_base = data_table.Grid_CO2_hrs_kg_(1:N)./907.185; % excess CO2 emissions in US ton
% ensure any time there is no excess load demand, 
% there is also no excess CO2 emissions
co2_t_base(ell_t_base == 0.0) = 0.0;

% total excess CO2 emissions for the baseline network
co2_emissions_baseline = sum(co2_t_base);

%% prepping storage arrays

cycles_naive = zeros(length(pack_sizes), length(hours_shift), table_len);
cycles_optimal = zeros(length(pack_sizes), length(hours_shift), table_len);
carbon_optim_naive = zeros(length(pack_sizes), length(hours_shift), table_len);
carbon_optim_optimal = zeros(length(pack_sizes), length(hours_shift), table_len);

%% pretty up output files
carbon_filename = ['./carbon_output' ...
                   '_mu_' num2str(mu, '%.2f') ...
                   '_P2Eratio_' num2str(power_to_energy_ratio, '%.2f') ...
                   '_SolarFarmSiteVariation_MedianTemp_outmat.mat']; % output file name
count_filename = ['./numCycle_output' ...
                  '_mu_' num2str(mu, '%.2f') ...
                  '_P2Eratio_' num2str(power_to_energy_ratio, '%.2f') ...
                  '_SolarFarmSiteVariation_MedianTemp_outmat.mat']; % output file name
variable_names = cell(13, 1);
variable_names(:) = {'CO2_intensity_t', ...
                     'Pdemand_total_MW', ...
                     'Ppv_total_MW', ...
                     'Pgrid_total_MW', ...
                     'Pgrid_to_load_MW', ...
                     'Pgrid_to_ESS_MW', ...
                     'Ppv_to_load_MW', ...
                     'Ppv_to_ESS_MW', ...
                     'Ppv_crtail_MW', ...
                     'Pess_chg_MW', ...
                     'Pess_dchg_MW', ...
                     'Eess_MWh', ...
                     'dchg_indicator'};

%% main loop

% optimization options
opts = optimoptions('intlinprog', ...
                    'Display', 'none', ...
                    'Heuristics', 'advanced');

% loop through all sites in region
for k=1:table_len

    pv_cap_t_base = table2array(pv_data_table(k, 6:end));
    pv_t = circshift(pv_cap_t_base(:), -5); % profiles are in UTC, shift to EST

    for s=1:length(hours_shift) % considering demand shifts
    
        hours_to_shift = hours_shift(s);
    
        ell_t = circshift(ell_t_base(:), hours_to_shift);
        co2_t = circshift(co2_t_base(:), hours_to_shift);

        nonzero_ell_t_idxs = find(ell_t ~= 0.0);

        zero_ell_t_idxs = find(ell_t == 0.0);
        zero_co2_t_idxs = find(co2_t == 0.0);
        zero_idxs = union(zero_ell_t_idxs, zero_co2_t_idxs);
    
        filename = ['./full_output' ...
                    '_shift_' num2str(hours_to_shift, '%.0f') ...
                    '_mu_' num2str(mu, '%.2f') ...
                    '_P2Eratio_' num2str(power_to_energy_ratio, '%.2f') ...
                    '_SolarFarmSite_' num2str(k, '%03.0f') ...
                    '_MedianTemp_outmat.xlsx']; % output file name
    
        for ps=1:length(pack_sizes)
    
            pack_size = pack_sizes(ps); % size of pack in MWh

            % maximum power in MW based on P2E ratio
            P_ESS = power_to_energy_ratio.*pack_size; 

            % ============= OPTIMAL PROBLEM =============
            
            % initialize the problem
            prob = optimproblem('ObjectiveSense', 'maximize');
    
            % (valid) trial solution - assume baseline network
            x0.g_t = ell_t;
            x0.g_load = ell_t;
            x0.g_ESS = zeros_vec;
            x0.pv_load = zeros_vec;
            x0.pv_ESS = zeros_vec;
            x0.pv_crtail = pv_t;
            x0.c_t = zeros_vec;
            x0.d_t = zeros_vec;
            x0.E_t = starting_SoE.*pack_size.*battery_present.*ones(N+1, 1);
            x0.dbin_t = ones(N, 1);
            
            % define decision variables
            g_t       = optimvar('g_t', N, 1, 'Type', 'continuous', 'LowerBound', 0.0);
            g_load    = optimvar('g_load', N, 1, 'Type', 'continuous', 'LowerBound', 0.0);
            g_ESS     = optimvar('g_ESS', N, 1, 'Type', 'continuous', 'LowerBound', 0.0);
            pv_load   = optimvar('pv_load', N, 1, 'Type', 'continuous', 'LowerBound', 0.0);
            pv_ESS    = optimvar('pv_ESS', N, 1, 'Type', 'continuous', 'LowerBound', 0.0);
            pv_crtail = optimvar('pv_crtail', N, 1, 'Type', 'continuous', 'LowerBound', 0.0);
            c_t       = optimvar('c_t', N, 1, 'Type', 'continuous', 'LowerBound', 0.0);
            d_t       = optimvar('d_t', N, 1, 'Type', 'continuous', 'LowerBound', 0.0);
            E_t       = optimvar('E_t', N+1, 1, 'Type', 'continuous', 'LowerBound', SOE_min.*pack_size.*battery_present, 'UpperBound', SOE_max.*pack_size.*battery_present);
            dbin_t    = optimvar('dbin_t', N, 1, 'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);
            
            % define objective function - total excess CO2 emissions
            % reduction
            prob.Objective = 1e2.*(1.0 - (...
                             sum(g_t(nonzero_ell_t_idxs).*co2_t(nonzero_ell_t_idxs)./ell_t(nonzero_ell_t_idxs)) ...
                             )./co2_emissions_baseline);
            
            % inequality constraints
            prob.Constraints.ineq1 = d_t <= dbin_t.*P_ESS.*battery_present;
            prob.Constraints.ineq2 = c_t <= (1 - dbin_t).*P_ESS.*battery_present;
            
            % equality constraints              
            prob.Constraints.eq1 = ell_t          == g_load + pv_load + d_t;
            prob.Constraints.eq2 = pv_t           == pv_load + pv_ESS + pv_crtail;
            prob.Constraints.eq3 = g_t            == g_ESS + g_load;
            prob.Constraints.eq4 = c_t            == pv_ESS + g_ESS;
            prob.Constraints.eq5 = E_t(2:end)     == E_t(1:end-1) + mu.*c_t - d_t;
            prob.Constraints.eq6 = E_t([1; end])  == starting_SoE.*pack_size.*battery_present;
            prob.Constraints.eq7 = g_t(zero_idxs) == 0.0;
            
            % solve the problem
            [sol,fval,exitflag,output] = solve(prob, x0, 'Solver', 'intlinprog', 'Options', opts);
    
            % write optimal solution to excel file
            output_array = [co2_t ell_t pv_t sol.g_t sol.g_load sol.g_ESS sol.pv_load sol.pv_ESS sol.pv_crtail sol.c_t sol.d_t sol.E_t(2:end) sol.dbin_t]; 
            
            output_table = array2table(output_array);
            output_table.Properties.VariableNames(1:13) = variable_names;
            
            writetable(output_table, filename, ...
                       'Sheet', [num2str(pack_size, '%.0f') ' MWh Battery - Optimal'], ...
                       'WriteVariableNames', true);

            % compute the number of cycles consumed in a year   
            cycles_optimal(ps, s, k) = compute_cycles([0.0; hours_in_year], sol.E_t./pack_size);

            % record optimal value
            carbon_optim_optimal(ps, s, k) = fval;
        
        end
    end
end

%% save cycles and optimal objective function values

save(count_filename, 'cycles_optimal');
save(carbon_filename, 'carbon_optim_optimal');

%% Helper functions

function hrs = return_hours_of_day(day_idxs)

    hrs = zeros(24, length(day_idxs));

    for i=1:length(day_idxs)

        n = day_idxs(i);

        hrs_low = 24.*(n-1) + 1;
        hrs_high = 24.*n;
        
        hrs(:,i) = (hrs_low:hrs_high)';
    end

end

function [num_cycles] = compute_cycles(time_idxs, SOE_t)
    
    dSOE_dt_full = gradient(SOE_t, time_idxs);

    SOE_chg = 0.0;
    SOE_dchg = 0.0;

    cycle_chg = 0;
    cycle_dchg = 0;

    for i = 1:length(dSOE_dt_full)

        dSOE_dt = dSOE_dt_full(i);

        if dSOE_dt > 0.0
            SOE_chg = SOE_chg + abs(dSOE_dt);
        else
            SOE_dchg = SOE_dchg + abs(dSOE_dt);
        end

        if SOE_chg >= 1.0
            cycle_chg = cycle_chg + 1;
            SOE_chg = SOE_chg - 1.0;
        end

        if SOE_dchg >= 1.0
            cycle_dchg = cycle_dchg + 1;
            SOE_dchg = SOE_dchg - 1.0;
        end

    end

    num_cycles = 0.5.*(cycle_chg + cycle_dchg);

end