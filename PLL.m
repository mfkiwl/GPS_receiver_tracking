function [f, phi] = PLL (I_sum, Q_sum, satellite_number, integrate_periods) % Costas loop
    phi = atan (Q_sum / I_sum);
    %fprintf (1, 'satellite_index = %i, phi = %i\n', satellite_number, phi);

    f = PLL_loop_filter (phi, satellite_number, integrate_periods);
end