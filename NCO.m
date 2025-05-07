function s = NCO (f, satellite_number, sample_rate)
    persistent phase_accumulator;
    if isempty (phase_accumulator)
        satellites_quantity = 31;
        phase_accumulator = zeros (satellites_quantity, 1);
    end

    s = exp (1j * phase_accumulator (satellite_number));
    
    phase_accumulator (satellite_number) = phase_accumulator (satellite_number) + 2 * pi * f / sample_rate;
    while phase_accumulator (satellite_number) > 2 * pi
        phase_accumulator (satellite_number) = phase_accumulator (satellite_number) - 2 * pi;
    end
end
