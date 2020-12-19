function [idx, C, indices, dists] = kcentergreedy(X, num_clusters)
    n = size(X, 1);
    indices = zeros(1, num_clusters);

    if num_clusters >= size(X, 1)
        idx = 1:size(X, 1);
        C = X;
        indices = 1:size(X, 1);
        dists = zeros(size(X, 1), 1);
        return
    end

    xsum = sum(X.^2, 2);

    indices(1) = randi(n, 1, 1);
    dists = pdist2_(X, X(indices(1), :), xsum);
    idx = ones(n, 1);

    for i = 2:num_clusters
        [~, new_id] = max(dists);
        dists_new = pdist2_(X, X(new_id, :), xsum);

        idx(dists > dists_new) = i;
        dists = min(dists, dists_new);
        indices(i) = new_id;
    end

    C = X(indices, :);
end

function out = pdist2_(x, y, xsum)
    ysum = sum(y.^2, 2);
    out = bsxfun(@minus, xsum, x * (2 * y)');
    out = bsxfun(@plus, ysum', out);
end
