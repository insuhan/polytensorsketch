function Z = get_rbf_features(X, degree, sketch_dim, method, sampling, num_clusters)

    if strcmpi(method, 'grr') || strcmpi(method, 'taylor') || strcmpi(method, 'cheby')
        func = @(x) exp(2 * x);
        xsum = sum(X.^2, 2);

        c_params.sampling = 'kcenter';
        c_params.num_clusters = 10;
        c_params.xmin = -max(xsum);
        c_params.xmax = max(xsum);
        c_params.xx = X;
        c_params.yy = X;
        c_params.m = sketch_dim;

        if strcmpi(method, 'taylor')
            c_params.anchor = mean(X(:));
            c_params.const = 2;
        end

        [coeffs, ~] = get_coeffs(func, degree, method, c_params);

        Z = poly_tensor_sketch(X, coeffs, degree, sketch_dim);
        Z = Z .* exp(-xsum);

    elseif strcmpi(method, 'rff')% Random Fourier Features for RBF kernel
        num_features = size(X, 2);
        num_samples = 1 + degree * sketch_dim;
        w = sqrt(2) * randn(num_samples, num_features);
        u = 2 * pi * rand(1, num_samples);
        Z = sqrt(2 / num_samples) * cos(X * w' + u);

    else
        error("Not a valid option.")
    end

end
