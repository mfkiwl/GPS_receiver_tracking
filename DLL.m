function shift = DLL (sample, LO_value, PRN_index, satellite_number, samples_per_chip)
    offset_samples = 3; % 3 / 5 = 0.6 chip

    persistent PRN I_E_FIFO I_P_FIFO I_L_FIFO Q_E_FIFO Q_P_FIFO Q_L_FIFO;
    if isempty (PRN) % first use of function
        satellites_quantity = 31;
        PRN = cacode (linspace (1, satellites_quantity, satellites_quantity), samples_per_chip); % create a local copy of PRN codes

        period_chips = 1023;
        period_samples = period_chips * samples_per_chip;
        I_E_FIFO = zeros (satellites_quantity, period_samples);
        I_P_FIFO = zeros (satellites_quantity, period_samples);
        I_L_FIFO = zeros (satellites_quantity, period_samples);
        Q_E_FIFO = zeros (satellites_quantity, period_samples);
        Q_P_FIFO = zeros (satellites_quantity, period_samples);
        Q_L_FIFO = zeros (satellites_quantity, period_samples);
    end
    
    %fprintf (1, 'satellite_index = %i, PRN_index = %i\n', satellite_number, PRN_index);

    if PRN_index + 3 > samples_per_chip
        PRN_index_E = PRN_index + 3 - samples_per_chip;
    else 
        PRN_index_E = PRN_index + 3;
    end

    if PRN_index - 3 < 1
        PRN_index_L = PRN_index - 3 + samples_per_chip;
    else 
        PRN_index_L = PRN_index - 3;
    end

    PRN_E = PRN (satellite_number, PRN_index_E);
    PRN_P = PRN (satellite_number, PRN_index);
    PRN_L = PRN (satellite_number, PRN_index_L);
    
    I_Q = sample * LO_value;
    I = real (I_Q);
    Q = imag (I_Q);

    I_E_new = I * PRN_E;
    I_P_new = I * PRN_P;
    I_L_new = I * PRN_L;

    Q_E_new = Q * PRN_E;
    Q_P_new = Q * PRN_P;
    Q_L_new = Q * PRN_L;

    % update FIFOs
    I_E_FIFO (satellite_number, 2 : end) = I_E_FIFO (satellite_number, 1 : end - 1);
    I_E_FIFO (satellite_number, 1) = I_E_new;

    I_P_FIFO (satellite_number, 2 : end) = I_P_FIFO (satellite_number, 1 : end - 1);
    I_P_FIFO (satellite_number, 1) = I_P_new;

    I_L_FIFO (satellite_number, 2 : end) = I_L_FIFO (satellite_number, 1 : end - 1);
    I_L_FIFO (satellite_number, 1) = I_L_new;

    Q_E_FIFO (satellite_number, 2 : end) = Q_E_FIFO (satellite_number, 1 : end - 1);
    Q_E_FIFO (satellite_number, 1) = Q_E_new;

    Q_P_FIFO (satellite_number, 2 : end) = Q_P_FIFO (satellite_number, 1 : end - 1);
    Q_P_FIFO (satellite_number, 1) = Q_P_new;

    Q_L_FIFO (satellite_number, 2 : end) = Q_L_FIFO (satellite_number, 1 : end - 1);
    Q_L_FIFO (satellite_number, 1) = Q_L_new;

    % sum values in FIFOs
    I_E = sum (I_E_FIFO (satellite_number, :));
    I_P = sum (I_P_FIFO (satellite_number, :));
    I_L = sum (I_L_FIFO (satellite_number, :));

    Q_E = sum (Q_E_FIFO (satellite_number, :));
    Q_P = sum (Q_P_FIFO (satellite_number, :));
    Q_L = sum (Q_L_FIFO (satellite_number, :));

    % discriminator
    if (I_E^2 + Q_E^2) + (I_L^2 + Q_L^2) > 0
        d = ((I_E^2 + Q_E^2) - (I_L^2 + Q_L^2)) / ((I_E^2 + Q_E^2) + (I_L^2 + Q_L^2));
    else
        d = 0;
        fprintf (1, 'satellite_index = %i, NaN avoided\n', satellite_number);
    end

    %fprintf (1, 'satellite_index = %i, d = %i\n', satellite_number, d);

    % shift calculation
    shift = round (d * 0.5 * samples_per_chip);
end