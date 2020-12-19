function TUp = tensor_sketch(U, m, p, inds_rand, sign_rand)
    k = 1;
    CU = fft(CountSketchMex(U, inds_rand(k, :), m, sign_rand(k, :), 1), [], 2);

    for k = 2:p
        CU = CU .* fft(CountSketchMex(U, inds_rand(k, :), m, sign_rand(k, :), 1), [], 2);
    end

    TUp = ifft(CU, [], 2);
end
