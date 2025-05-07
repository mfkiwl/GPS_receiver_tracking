function f = PLL (I_sum, Q_sum, satellite_number) % Costas loop
    phi = atan (Q_sum / I_sum);

    f = loop_filter (phi, satellite_number);
end