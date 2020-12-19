function [c, params] = get_coeffs(func, degree, method, params)

    if strcmp(method, 'taylor')
        [c, params] = taylor_coeffs(func, degree, method, params);

    elseif strcmp(method, 'cheby')
        xmin = params.xmin; xmax = params.xmax;
        g = @(x) x .* ((xmax - xmin) / 2) + ((xmax + xmin) / 2);
        xk = cos(pi * ((0:degree)' + 0.5) / (degree + 1)); % zeros of Tn(x)
        fk = func(g(xk)); % target function, [-1,1]->[del,1-del]
        Tk = ones(degree + 1, degree + 1); Tk(:, 2) = xk; % init recursion

        for i = 2:degree
            Tk(:, i + 1) = 2 * xk .* Tk(:, i) - Tk(:, i - 1); % evaluate polynomials
        end

        c = 2 / (degree + 1) * (Tk' * fk); c(1) = c(1) / 2; % compute Chebyshev coefficients
        params.cheby_coeffs = c;
        % Chebyshev -> usual coefficients
        pp = zeros(degree + 1, degree + 1);
        pp(:, 1) = [1; zeros(degree, 1)];
        pp(:, 2) = [0; 1; zeros(degree - 1, 1)];

        for j = 2:degree
            pp(:, j + 1) = [0; 2 * pp(1:degree, j)] - pp(:, j - 1);
        end

        % p(ginv(x)) -> monomial based coefficients
        a = (2 / (xmax - xmin)); b = -(xmax + xmin) / (xmax - xmin);
        pp2 = zeros(degree + 1, degree + 1);
        pp2(1, :) = b.^(0:degree);

        for i = 1:degree
            tmp = 1;

            for j = 0:i - 1
                tmp = tmp .* ((i:degree) - j);
            end

            pp2(i + 1, i + 1:end) = tmp .* (b.^(0:degree - i)) * a.^i / factorial(i);
        end

        c = (pp2 * pp * reshape(c, [], 1))';

    elseif strcmpi(method, 'grr') % coefficients based-on generalized ridge regression 
        W = get_W(params, degree);
        xmin = params.xmin;
        xmax = params.xmax;
        g = @(x) x .* (2 / (xmax - xmin)) - (xmax + xmin) / (xmax - xmin);

        if ~isfield(params, 'sampling')
            params.sampling = 'kcenter';
            params.num_clusters = 10;
        end

        if strcmp(params.sampling, 'uniform')
            [x, rescale] = get_entries_uniform(params);
        elseif contains(params.sampling, {'kmedoids', 'kmeans', 'kcenter', 'uniform'})
            [x, rescale] = get_entries_coreset(params);
        elseif contains(params.sampling, 'exact')
            x = reshape(params.xx * params.xx', [], 1);
            rescale = ones(size(x));
        end

        gx = g(x);
        Xc = ones(length(x), degree + 1);
        Xc(:, 2) = gx;

        for k = 2:degree
            Xc(:, k + 1) = 2.0 * gx .* Xc(:, k) - Xc(:, k - 1);
        end

        a = (2.0 / (xmax - xmin));
        b = -(xmax + xmin) / (xmax - xmin);
        P = get_transform(degree, a, b);
        WP = W * P;
        f = func(x);

        XXWW = (Xc' * (rescale .* Xc)) + WP' * WP;
        Xf = Xc' * (rescale .* f);
        c_ = XXWW \ Xf;
        c = P * c_;

    elseif strcmp(method, 'rff') || strcmp(method, 'nystrom')
        c = zeros(1, degree + 1);
    end

end

function W = get_W(params, degree)

    if norm(params.xx - params.yy, 'fro') < 1e-5
        U = params.xx;
        m = params.m;
        w = zeros(1, degree + 1);

        for i = 2:degree + 1
            w(i) = sqrt(degree * (2 + 3^i) / m) * sum(sum(U.^2, 2).^i);
        end

        W = diag(w);
        return
    end

    U = params.xx;
    V = params.yy;
    m = params.m;
    W = zeros(degree + 1, degree + 1);

    for i = 2:degree + 1
        W(i, i) = sqrt(degree * (2 + 3^i) * sum(sum(U.^2, 2).^i) * sum(sum(V.^2, 2).^i) / m);
    end

end

function p = get_transform(degree, a, b)
    p1 = zeros(degree + 1, degree + 1);
    p1(:, 1) = [1; zeros(degree, 1)];
    p1(:, 2) = [0; 1; zeros(degree - 1, 1)];

    for j = 2:degree
        p1(:, j + 1) = [0; 2 * p1(1:degree, j)] - p1(:, j - 1);
    end

    p = zeros(degree + 1, degree + 1);
    p(1, :) = b.^(0:degree);

    for i = 1:degree
        tmp = 1;

        for j = 0:i - 1
            tmp = tmp .* ((i:degree) - j);
        end

        p(i + 1, i + 1:end) = tmp .* (b.^(0:degree - i)) * a.^i / factorial(i);
    end

    p = p * p1;
end
