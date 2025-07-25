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