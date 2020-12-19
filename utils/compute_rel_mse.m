function err = compute_rel_mse(X, X_true)
    err = norm( (X - X_true) ./ X_true, 'fro').^2 / numel(X_true);
end

