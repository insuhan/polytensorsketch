function X = load_dataset(dataset, num_points, num_features)

    if ismember(dataset, {'satimage', 'segment'})
        A = load(sprintf('data/%s.mat', dataset));
        X = A.feats;
    elseif ismember(dataset, 'synthetic')
        X = randn(num_points, num_features);
    else
        error(sprintf('%s is not supported yet.', dataset));
    end

end
