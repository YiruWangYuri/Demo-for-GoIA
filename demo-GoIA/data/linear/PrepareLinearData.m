function [data, para] = PrepareLinearData(N, outlierRatio, dim, seed, noise_in, noise_out)
    FlagofBlance = 0;
    [X_Data, theta, th_guide] = genSyntheticData(N, noise_in, noise_out, outlierRatio, dim, FlagofBlance, seed);
    data = X_Data(:, 1:end-1)';
    para = [N, dim, noise_in, outlierRatio, th_guide, theta'];    
end