function f = loop_filter (phi_new, satellite_number) % Second-order phase lock loop filter
    persistent memory C_1 C_2;
    if isempty (memory)
        % Execute this code for a single time to calculate filter
        % coefficients C_1 and C_2
        % K_0 * K_d = 1
        zeta = 20; % damping ratio
        B_L = 10000; % noise bandwidth in Hz
        periods_per_second = 1000;
        T = 1 / periods_per_second; % loop filter sampling time in s
        omega_n = 8 * zeta * B_L / (4 * zeta^2 + 1); % loop filter natural frequency 
        
        C_1 = 8 * zeta * omega_n * T / (4 + 4 * zeta * omega_n * T + (omega_n * T)^2);
        C_2 = 4 * (omega_n * T)^2 / (4 + 4 * zeta * omega_n * T + (omega_n * T)^2);
        fprintf (1, 'C_1 = %d, C_2 = %d\n', C_1, C_2);

        satellites_quantity = 31;
        memory = zeros (satellites_quantity, 1);
    end

    branch_1 = phi_new * C_1;
    branch_2 = phi_new * C_2 + memory (satellite_number);
    memory (satellite_number) = branch_2;
    f = branch_1 + branch_2;
end