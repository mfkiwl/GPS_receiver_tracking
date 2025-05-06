function f = PLL (sample, LO_value, code_value, satellite_number) % Costas loop
    I_Q = sample * code_value * LO_value;
    I = real (I_Q);
    Q = imag (I_Q);

    [I_LPF, Q_LPF] = LPF (I, Q, satellite_number);

    phi = atan (Q_LPF / I_LPF);

    f = loop_filter (phi, satellite_number);
end