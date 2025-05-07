function shift = DLL (I_E_sum, Q_E_sum, I_L_sum, Q_L_sum, satellite_number, samples_per_chip)
    %fprintf (1, 'satellite_index = %i, PRN_index = %i\n', satellite_number, PRN_index);

    % discriminator
    if (I_E_sum^2 + Q_E_sum^2) + (I_L_sum^2 + Q_L_sum^2) > 0
        d = ((I_E_sum^2 + Q_E_sum^2) - (I_L_sum^2 + Q_L_sum^2)) / ((I_E_sum^2 + Q_E_sum^2) + (I_L_sum^2 + Q_L_sum^2));
    else
        d = 0;
        fprintf (1, 'satellite_index = %i, NaN avoided\n', satellite_number);
    end

    %fprintf (1, 'satellite_index = %i, d = %i\n', satellite_number, d);

    % shift calculation
    shift = round (d * 0.5 * samples_per_chip);
end