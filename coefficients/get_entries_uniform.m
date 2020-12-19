function [vals, rescale] = get_entries_uniform(params)
    ns = params.ns;
    X = params.xx;
    n = size(X, 1);

    idx = randsample(n * n, ns);

    [xidx, yidx] = ind2sub([n, n], idx);
    vals = zeros(ns, 1);

    for k = 1:length(xidx)
        vals(k) = X(xidx(k), :) * X(yidx(k), :)';
    end

    rescale = ones(ns, 1) * (n^2 / ns);
end
