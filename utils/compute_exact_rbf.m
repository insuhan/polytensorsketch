function Z = compute_exact_rbf(X)
    exp_l2norm = exp(-sum(X.^2, 2));
    exp_XXT = exp(2.0 * (X * X'));
    Z = exp_l2norm .* (exp_l2norm .* exp_XXT)';
end
