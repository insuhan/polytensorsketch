function Z = poly_tensor_sketch(X, coeffs, degree, sketch_dim)
    [n, d] = size(X);
    out_dim = 1 + degree * sketch_dim;
    Z = zeros(n, out_dim);

    c_sqrt = sqrt(coeffs);
    Z(:, 1) = c_sqrt(1);
    last_idx = 1;

    for k = 1:degree
        begin_idx = last_idx + 1;
        last_idx = begin_idx + sketch_dim - 1;

        rand_hash = int64(randi(sketch_dim, k, d));
        rand_sign = sign(randn(k, d));

        feats = tensor_sketch(X, sketch_dim, k, rand_hash, rand_sign);
        Z(:, begin_idx:last_idx) = c_sqrt(k + 1) * feats;
    end

end
