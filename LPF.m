function [I_output, Q_output] = LPF (I_new, Q_new, satellite_number)
    persistent I_array Q_array h;
    if isempty (h) % first use of function
        cutoff_frequency = 700; %Hz
        sample_rate = 5 * 10^6 * 1.023; %Hz
        filter_order = 6000;
        h = fir1 (filter_order, cutoff_frequency / (0.5 * sample_rate), "low");
    
        satellites_quantity = 31;
        I_array = zeros (satellites_quantity, length (h));
        Q_array = zeros (satellites_quantity, length (h));
    end
    
    % update FIFOs
    I_array (satellite_number, 2 : end) = I_array (satellite_number, 1 : end - 1);
    I_array (satellite_number, 1) = I_new;

    Q_array (satellite_number, 2 : end) = Q_array (satellite_number, 1 : end - 1);
    Q_array (satellite_number, 1) = Q_new;
    % apply FIR LPF
    I_output = FIR_filter (I_array (satellite_number, :), h);
    Q_output = FIR_filter (Q_array (satellite_number, :), h);
end