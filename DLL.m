function shift = DLL (I_E_sum, Q_E_sum, I_L_sum, Q_L_sum, satellite_index, integrate_periods)
    % discriminator
    if (I_E_sum^2 + Q_E_sum^2) + (I_L_sum^2 + Q_L_sum^2) > 0
        D = ((I_E_sum^2 + Q_E_sum^2) - (I_L_sum^2 + Q_L_sum^2)) / ((I_E_sum^2 + Q_E_sum^2) + (I_L_sum^2 + Q_L_sum^2));
    else
        D = 0;
        fprintf (1, 'satellite_index = %i, NaN avoided\n', satellite_index);
    end

    %fprintf (1, 'satellite_index = %i, d = %i\n', satellite_number, d);

    % shift calculation
    %shift = round (D * 0.5 * samples_per_chip);
    shift = round (DLL_loop_filter (D, satellite_index, integrate_periods));
end