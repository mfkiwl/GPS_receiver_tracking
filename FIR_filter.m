function output = FIR_filter (x, h)
    sum = 0;
    size = length (h);

    for k = 1 : 1 : size
        sum = sum + h (k) * x (size + 1 - k);
    end

    output = sum;
end