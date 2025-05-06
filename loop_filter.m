function f = loop_filter (phi_new, satellite_number) % Second-order phase lock loop filter
    % K_0 * K_d = 1
    zeta = 0.7; % damping ratio
    B_L = 30; % noise bandwidth in Hz
    T = 1 / (5 * 10^6 * 1.023); % loop filter sampling time
    omega_n = 8 * zeta * B_L / (4 * zeta^2 + 1); % loop filter natural frequency 

    persistent memory C_1 C_2;
    if isempty (memory)
        satellites_quantity = 31;
        
        C_1 = 8 * zeta * omega_n * T / (4 + 4 * zeta * omega_n * T + (omega_n * T)^2);
        C_2 = 4 * (omega_n * T)^2 / (4 + 4 * zeta * omega_n * T + (omega_n * T)^2);
        fprintf (1, 'C_1 = %d, C_2 = %d\n', C_1, C_2);

        memory = zeros (satellites_quantity, 1);
    end

    branch_1 = phi_new * C_1;
    branch_2 = phi_new * C_2 + memory (satellite_number);
    memory (satellite_number) = branch_2;
    f = branch_1 + branch_2;
end