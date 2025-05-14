function shift = DLL_loop_filter (D_new, satellite_number, integrate_periods) % Second-order phase lock loop filter
    persistent memory C_1 C_2;
    if isempty (memory)
        % Execute this code for a single time to calculate filter
        % coefficients C_1 and C_2
        % K_0 * K_d = 1
        zeta = 1; % damping ratio
        B_L = 2; % noise bandwidth in Hz

        periods_per_second = 1000;
        %sample_rate = 5115000;
        T = integrate_periods / periods_per_second; % loop filter sampling time in s
        %T = 1 / sample_rate;

        omega_n = 8 * zeta * B_L / (4 * zeta^2 + 1); % loop filter natural frequency 
        
        %C_1 = 8 * zeta * omega_n * T / (4 + 4 * zeta * omega_n * T + (omega_n * T)^2);
        %C_2 = 4 * (omega_n * T)^2 / (4 + 4 * zeta * omega_n * T + (omega_n * T)^2);

        C_1 = (2 * zeta / omega_n) / (1 / (omega_n^2));
        C_2 = T / (1 / (omega_n^2));
        fprintf (1, 'DLL: C_1 = %d, C_2 = %d\n', C_1, C_2);

        satellites_quantity = 31;
        memory = zeros (satellites_quantity, 1);
    end

    branch_1 = D_new * C_1;
    branch_2 = D_new * C_2 + memory (satellite_number);
    memory (satellite_number) = branch_2;
    shift = branch_1 + branch_2;
end