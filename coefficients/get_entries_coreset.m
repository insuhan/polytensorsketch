function [entries, rescale] = get_entries_coreset(params)
    assert(contains(params.sampling, {'kmedoids', 'kmeans', 'kcenter', 'uniform'}))

    if ~isfield(params, 'num_clusters')
        fprintf("default number of clusters: %d\n", 10);
        num_clusters = 10;
    else
        num_clusters = params.num_clusters;
    end

    U = params.xx;
    V = params.yy;

    if norm(U - V, 'fro') < 1e-10 % symmetric case such that U = V

        if strcmp(params.sampling, 'kmedoids')
            [idx_U, U_C] = kmedoids(U, num_clusters, 'distance', 'euclidean');
        elseif strcmp(params.sampling, 'kmeans')
            [idx_U, U_C] = kmeans(U, num_clusters);
        elseif strcmp(params.sampling, 'kcenter')
            [idx_U, U_C, ~] = kcenter_greedy(U, num_clusters);
        elseif strcmp(params.sampling, 'uniform')
            params.ns = num_clusters * size(U, 1);
            [entries, rescale] = get_entries_uniform(params);
            return
        end

        entries = reshape(U_C * V', [], 1);
        rescale = reshape(histcounts(idx_U, (1:num_clusters + 1) - 0.5)' * ones(1, size(V, 1)), [], 1);
    else % nonsymmetric case
        [idx_U, U_C, U_dist_sum] = kmedoids(U, num_clusters, 'distance', 'euclidean', 'Algorithm', 'large');
        [idx_V, V_C, V_dist_sum] = kmedoids(V, num_clusters, 'distance', 'euclidean', 'Algorithm', 'large');

        U_l2norm_sum = sqrt(sum(U.^2, 2));
        V_l2norm_sum = sqrt(sum(V.^2, 2));

        if (U_l2norm_sum / V_l2norm_sum) < sum(U_dist_sum) / sum(V_dist_sum)
            entries = reshape(U_C * V', [], 1);
            rescale = reshape(histcounts(idx_U)' * ones(1, size(V, 1)), [], 1);
        else
            entries = reshape(U * V_C', [], 1);
            rescale = reshape(ones(size(U, 1), 1) * histcounts(idx_V), [], 1);
        end

    end

end
