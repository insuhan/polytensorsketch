function run_kernel_approximation(dataset)

    if nargin < 1
        dataset = 'synthetic';
        X = load_dataset(dataset, 1000, 50);
    elseif strcmpi(dataset, 'synthetic')
        X = load_dataset(dataset, 1000, 50);
    else
        X = load_dataset(dataset, -1, -1);
    end

    fprintf("%s dataset is loaded, n=%d, d=%d\n", dataset, size(X,1), size(X,2));
    sc = sqrt(10);

    % Normalize the input matrix
    [n, d] = size(X);
    X = X / sqrt(d) / sc;

    degree = 10;
    sketch_dim = 10;

    % Compute the exact RBF kernel matrix
    K_exact = compute_exact_rbf(X);

    % Approximate RBF kernel using polynomial tensor sketch with coreset coefficients
    method = 'grr';
    sampling = 'kcenter';
    num_clusters = 10;
    Z_coreset = get_rbf_features(X, degree, sketch_dim, method, sampling, num_clusters);
    K_coreset = Z_coreset * Z_coreset';
    err_coreset = compute_rel_mse(K_coreset, K_exact);
    fprintf("PTS (coreset) error: %.6f\n", err_coreset);

    % Approximate RBF kernel via Random Fourier Features
    Z_rff = get_rbf_features(X, degree, sketch_dim, 'rff', -1, -1);
    K_rff = Z_rff * Z_rff';
    err_rff = compute_rel_mse(K_rff, K_exact);
    fprintf("RFF           error: %.6f\n", err_rff);

end