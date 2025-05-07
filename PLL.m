function [f, phi] = PLL (I_sum, Q_sum, satellite_number) % Costas loop
    phi = atan (Q_sum / I_sum);
    %fprintf (1, 'satellite_index = %i, phi = %i\n', satellite_number, phi);

    f = loop_filter (phi, satellite_number);
end