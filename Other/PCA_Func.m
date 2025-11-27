function [Accy_PCA, best_matches] = PCA_Func(x, y, img_size, train_files, test_files, save_folder, true_labels)
    x = x'; y = y';
    [a, b] = size(x);
    CovData = cov(x);
    [U, ~, ~] = svd(CovData);
    maxComp = min(50, b);
    Accy_PCA = zeros(1, maxComp); % preallocate
    best_matches = zeros(a, 1); % initialize

    for n = 1:maxComp
        U_n = U(:, 1:n);
        PCA_train = x * U_n;
        PCA_test  = y * U_n;

        NN_PCA = zeros(a, 1);
        for i = 1:a
            dist_PCA = vecnorm(PCA_test(i,:) - PCA_train, 2, 2);
            [~, idx] = min(dist_PCA);
            NN_PCA(i) = idx;
            if n == 10
                best_matches(i) = idx;
            end
        end

        predicted_labels = true_labels(NN_PCA)';
        Accy_PCA(n) = sum(predicted_labels == true_labels) / a;
    end

    figure();
    plot(1:maxComp, Accy_PCA, 'ko-', 'LineWidth', 2);
    xlabel('Number of Principal Components');
    ylabel('Classification Accuracy');
    title('PCA + Nearest Neighbor Accuracy');

    % Save images for best matches when n = 10
    for k = 1:a
        test_img = imread(test_files{k});
        imwrite(imresize(test_img, img_size), ...
            fullfile(save_folder, sprintf('Test_%04d.tif', k)));

        match_img = imread(train_files{best_matches(k)});
        imwrite(imresize(match_img, img_size), ...
            fullfile(save_folder, sprintf('Matched_%04d.tif', k)));
    end
end
