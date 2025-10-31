% =========================================================================
%
% ResNet Model Training Pipeline (PNIPAM-DNA-DiYO-1)
%
% Overview:
%   This script trains a ResNet-18 deep learning model for the
%   classification of single-cell images containing PNIPAM-DNA-DiYO-1
%   microparticles. The model classifies cells into three categories:
%   single-particle, multi-particle, or no-particle, based on features
%   extracted from segmented cell images.
%
% Input Files:
%   - FITC channel images (8-bit TIFF, format: ParticleType_##.tif)
%   - DAPI channel images (8-bit TIFF, format: ParticleType_##.tif)
%   - Brightfield (phase contrast) images (8-bit TIFF, format: ParticleType_##.tif)
%   *Note: All images must be pre-processed and split into separate
%          channels using ImageJ/Fiji before training.
%
% Segmentation and Masking:
%   Binary segmentation masks are generated from:
%       - Brightfield: cell boundaries
%       - DAPI: nucleus location
%       - FITC: particle signals
%   These masks are used to define ROIs and prepare labeled training data
%   for classification.
%
% Model Training (3-Class Classification):
%   - The ResNet-18 model is trained on segmented cell images containing
%     PNIPAM-DNA-DiYO-1 particles.
%   - Each image is assigned to one of three classes:
%       1) Single-particle
%       2) Multi-particle
%       3) No-particle
%   - Data augmentation (rotation, scaling, reflection) is applied to
%     enhance model robustness and generalization.
%
% Output Files:
%   1) Labeled cell images organized into three categories
%      (Single-particle, Multi-particle, No-particle)
%   2) Trained ResNet-18 classification model (.mat file)
%   3) Validation results, including accuracy/loss curves and confusion matrix
%
% Requirements:
%   - MATLAB R2023a or later
%   - Image Processing Toolbox
%   - Deep Learning Toolbox
%   - Pre-trained ResNet-18 (Deep Learning Toolbox Model for ResNet-18 Network)
%   - ImageJ/Fiji (for initial channel decomposition)
%
% =========================================================================

clc; clear; warning('off','all'); close all

% Directory (Example directory & Pathes)
BF_path    = 'D:/Image Process DNA/Training (ResNet)/PNI + DNA/BF';    mkdir(BF_path);
Green_path = 'D:/Image Process DNA/Training (ResNet)/PNI + DNA/Green'; mkdir(Green_path);

BFMask_path    = 'D:/Image Process DNA/Training (ResNet)/PNI + DNA/BF Mask';    mkdir(BFMask_path);
GreenMask_path = 'D:/Image Process DNA/Training (ResNet)/PNI + DNA/Green Mask'; mkdir(GreenMask_path);

CellBFSingle_path    = 'D:/Image Process DNA/Training (ResNet)/PNI + DNA/Single Particle/BF';      mkdir(CellBFSingle_path);
CellGreenSingle_path = 'D:/Image Process DNA/Training (ResNet)/PNI + DNA/Single Particle/Green';   mkdir(CellGreenSingle_path);
CellBFMulti_path     = 'D:/Image Process DNA/Training (ResNet)/PNI + DNA/Multiple Particle/BF';    mkdir(CellBFMulti_path);
CellGreenMulti_path  = 'D:/Image Process DNA/Training (ResNet)/PNI + DNA/Multiple Particle/Green'; mkdir(CellGreenMulti_path);
CellBFNone_path      = 'D:/Image Process DNA/Training (ResNet)/PNI + DNA/No Particle/BF';          mkdir(CellBFNone_path);
CellGreenNone_path   = 'D:/Image Process DNA/Training (ResNet)/PNI + DNA/No Particle/Green';       mkdir(CellGreenNone_path);

base_path = 'D:/Image Process DNA/Conv Watershed/PNI + DNA';

ImgCellSingle_count = 0;
ImgCellMulti_count  = 0;
ImgCellNone_count   = 0;

% ========================================================================== 
% 
%                   Cell segmentation
%
% ==========================================================================

phase_dir = fullfile(base_path, 'Phase');
blue_dir  = fullfile(base_path, 'Blue');
green_dir  = fullfile(base_path, 'Green');

% Image file sorting
phase_files = dir(fullfile(phase_dir, '*.tif'));
names = {phase_files.name};
nums = regexp(names, '\d+', 'match');
nums = cellfun(@(x) str2double(x{end}), nums);    
[~, idx] = sort(nums);
phase_files = phase_files(idx);

% Obtaining contrast & noise parameter
contrast_list = [];
noise_list = [];

for a = 1:length(phase_files)
    img = im2double(imread(fullfile(phase_dir, phase_files(a).name)));
    contrast_list(end+1) = std2(img);
    noise_list(end+1) = mean2(stdfilt(img));
end

% Threshold for Processing Images
Contrast_threshold = prctile(contrast_list, 20);
Noise_threshold = prctile(noise_list, 75);

for b = 1:length(phase_files)
% Image Acquisition
    phase_filename = phase_files(b).name;
    blue_filename  = strrep(phase_filename, 'Phase', 'Blue');
    green_filename  = strrep(phase_filename, 'Phase', 'Green');
    Img_1 = fullfile(phase_dir, phase_filename);
    Img_2 = fullfile(blue_dir, blue_filename);
    Img_3 = fullfile(green_dir, green_filename);

    if ~isfile(Img_1) || ~isfile(Img_2) || ~isfile(Img_3)
        continue;
    end

% Image Data Aquisition
    BF = imread(Img_1);     % Bright field image
    Blue = imread(Img_2);   % DAPI signal image
    Green = imread(Img_3);  % FITC signal image

% Parameters (Convolution)

    % For Bright-field
    CF1 = fspecial('sobel'); Trans_CF1 = CF1';
    Conv_Threshold_BF = 0.3;

    % For Nucleus (DAPI signal images)
    CF2 = fspecial('gaussian'); Trans_CF2 = CF2';
    Conv_Threshold_Nuc = 1.1;

    % Threshold for pre-processing images
    BW_noise_threshold = 0.15;

% Image Processing (Convolution Mask)

    % Bright field image mask
    I_BF = im2double(BF);
    img_contrast = std2(I_BF);
    img_noise = mean(mean(stdfilt(I_BF)));

    if img_noise > Noise_threshold
        sigma = min(2.0, 5 * (img_noise - Noise_threshold));
        I_BF = imgaussfilt(I_BF, sigma);
    end

    if img_contrast < Contrast_threshold
        stretch_range = max(0.01, 0.1 * (1 - img_contrast / Contrast_threshold));
        I_BF = imadjust(I_BF, stretchlim(I_BF, [stretch_range 1 - stretch_range]), [0 1]);
    end 

    I_med_BF = medfilt2(I_BF, [5 5]);
    I_sm = imbilatfilt(I_med_BF,std2(BF)^2, 0.5);

    edges_BF = sqrt((conv2(I_sm, CF1, 'same')).^2 + (conv2(I_sm, Trans_CF1, 'same')).^2);
    T_BF     = graythresh(edges_BF);
    bw_BF    = imbinarize(edges_BF, T_BF * Conv_Threshold_BF);
    bw_noise = mean(mean(stdfilt(bw_BF)));

    if bw_noise > BW_noise_threshold
        bw_BF = bwareaopen(bw_BF, 50);
    end

    bw_BF_clean = imclose(bw_BF,strel('disk',5));
    bw_BF_clean = imopen(bw_BF_clean, ones(5,5));
    bw_BF_clean = imfill(bw_BF_clean,'holes');
    smoothed = imgaussfilt(im2double(bw_BF_clean), 6);
    smoothed = smoothed > 0.5;
    bw_BF_clean = imerode(smoothed,strel('disk',5));
    mask_clean = bwareaopen(bw_BF_clean, 10);
    mask_clean = imclose(mask_clean, strel('disk', 2));
    mask_blur = imgaussfilt(im2double(mask_clean), 3);
    mask_smoothed = mask_blur > 0.7;
    mask_clean = imerode(mask_smoothed, strel('disk', 1));

    % Nucleus image mask
    I_Nuc = im2double(Blue);
    I_med_Nuc = imgaussfilt(I_Nuc,0.8);
    edges_Nuc = sqrt((conv2(I_med_Nuc, CF2, 'same')).^2 + (conv2(I_med_Nuc, Trans_CF2, 'same')).^2);
    T_Nuc = graythresh(edges_Nuc);
    bw_Nuc = imbinarize(edges_Nuc, T_Nuc * Conv_Threshold_Nuc);
    bw_Nuc_clean = imclose(bw_Nuc,strel('disk',10));
    bw_Nuc_clean = imdilate(bw_Nuc_clean, strel('disk', 8));
    bw_Nuc_clean = imopen(bw_Nuc_clean, ones(5,5));
    bw_Nuc_clean = imfill(bw_Nuc_clean,'holes');
    bw_Nuc_clean = imerode(bw_Nuc_clean,strel('diamond',7));

    % Microparticle image mask
    I_green = imbinarize(Green);
    D = bwdist(~I_green);

    bw_Green_clean = imextendedmax(D,2);
    bw_Green_clean = imdilate(bw_Green_clean, strel('disk',10));
    bw_Green_clean = imclearborder(bw_Green_clean);
    bw_Par = imerode(bw_Green_clean, strel('disk', 5));

    label_Mask = bwlabel(mask_clean);
    label_Nuc = bwlabel(bw_Nuc_clean);
    n_cells = max(label_Mask(:));

% Visual inspection: segmented figure (Overlay bright field w/ mask image)

    label_mask = bwlabel(bw_BF_clean);
    cent_num = regionprops(label_mask,'Centroid');

    for bb = 1:numel(cent_num)
        text(cent_num(bb).Centroid(1),cent_num(bb).Centroid(2),num2str(bb),'color','r')
    end

    figure()
    imshow(imoverlay(BF,bwperim(bw_BF_clean)))
    figure();
    imshow(imoverlay(imoverlay(imoverlay(mat2gray(BF),mask_clean,'blue'),bw_Nuc_clean,'red'),bw_Par,'green'));
    pause(1); close;

% ========================================================================== 
% 
%                   Separation of ROI
%
% ==========================================================================

% ROI identification / isolation
    for j = 1:n_cells
        cell_mask = (label_Mask == j);
        overlapped = label_Nuc .* cell_mask;
        unique_labels = unique(overlapped);
        unique_labels(unique_labels == 0) = [];
        nuclei_count = numel(unique_labels);

        [row1, col1] = find(cell_mask);
        len = max(row1) - min(row1) + 1;
        breadth = max(col1) - min(col1) + 1;
        target1 = uint8(zeros([len breadth]));
        target2 = uint8(zeros([len breadth]));
        target3 = zeros([len breadth]);
        target4 = zeros([len breadth]);
        target5 = zeros([len breadth]);
        sy = min(col1) - 1;
        sx = min(row1) - 1;


            for k = 1:length(row1)
                x = row1(k) - sx;
                y = col1(k) - sy;
                target1(x, y) = BF(row1(k), col1(k));
                target2(x, y) = Green(row1(k), col1(k));
                target3(x, y) = mask_clean(row1(k), col1(k));
                target4(x, y) = bw_Nuc_clean(row1(k), col1(k));
                target5(x, y) = bw_Par(row1(k), col1(k));
            end

        label_GreenMask_single = bwlabel(target5);
        num_particles = max(label_GreenMask_single(:));

    % Isolatung single cell (based on the number of nucleus in ROI)
    % Storing in classified directory
        if nuclei_count == 1
            nuc_erode = imerode(target4,strel('diamond',3));

            if num_particles == 1 
                imwrite(target1, fullfile(CellBFSingle_path, sprintf('%d.tif', ImgCellSingle_count + 1)));
                imwrite(target2, fullfile(CellGreenSingle_path, sprintf('%d.tif', ImgCellSingle_count + 1))); % Optional
                ImgCellSingle_count = ImgCellSingle_count + 1;
            elseif num_particles > 1
                imwrite(target1, fullfile(CellBFMulti_path, sprintf('%d.tif', ImgCellMulti_count + 1)));
                imwrite(target2, fullfile(CellGreenMulti_path, sprintf('%d.tif', ImgCellMulti_count + 1))); % Optional
                ImgCellMulti_count = ImgCellMulti_count + 1;
            else
                imwrite(target1, fullfile(CellBFNone_path, sprintf('%d.tif', ImgCellNone_count + 1)));
                imwrite(target2, fullfile(CellGreenNone_path, sprintf('%d.tif', ImgCellNone_count + 1))); % Optional
                ImgCellNone_count = ImgCellNone_count + 1;
            end

    % Segmenting aggregated sells
        elseif nuclei_count > 1
            I = im2double(target1);
            I = adapthisteq(I, 'ClipLimit', 0.05, 'Distribution', 'exponential', ...
                            'Alpha', 10, 'NBins', 200);
            I = imsharpen(I, 'Radius', 1.5, 'Amount', 2.5, 'Threshold', 0.01);
            I = imgaussfilt(I, 1.8);

            G = imgradient(I, 'sobel');
            G = imbilatfilt(G,std2(G)^2, 3);          
            G = imadjust(G, stretchlim(G), []);
            G = imcomplement(G);           

            mask_crop = imerode(target3 > 0, strel('disk', 3)); 
            nuclei    = imerode(target4 > 0, strel('disk', 3)); 

            G2 = imimposemin(G, nuclei);
            L  = watershed(G2);
            L(~mask_crop) = 0;
            L  = bwlabel(L > 0); 

            L_refined = zeros(size(L));

            for l = 1:max(L(:))
                region = (L == l);
                if nnz(region) < 100
                    continue;
                end
                initMask = imerode(region, strel('disk', 2));
                G_smooth = imgaussfilt(G, 5);
                refined_mask = activecontour(G_smooth, initMask, 8, 'edge');
                refined_mask_dilate = imdilate(refined_mask, strel('disk', 2));
                L_refined(refined_mask_dilate) = l;
            end

            label_refined = bwlabel(L_refined);

        % Visual validation (Overlay bright field w/ watershed image)
            % figure();
            % subplot(2,3,1); imshow(target1); title('Original BF Patch');
            % subplot(2,3,2); imshow(label2rgb(L, 'jet', 'k', 'shuffle')); title('Watershed Segmentation');
            % subplot(2,3,3); imshow(imoverlay(mat2gray(target1), bwperim(L_refined > 0), 'cyan')); title('Refined Edges (Active Contour)');
            % subplot(2,3,4); imshow(imoverlay(mat2gray(target1), target3, 'magenta')); title('Overlay with Nuclei');
            % subplot(2,3,5); imshow(imoverlay(mat2gray(target1), target4, 'magenta')); title('Overlay with Mask');
            % subplot(2,3,6); imshow(imoverlay(mat2gray(target1), target2, 'green')); title('Overlay with Partice');

        % Isolatung aggregated cell (based on the number of nucleus in ROI)
        % Storing in classified directory
            for m = 1:max(label_refined(:))
                [row2, col2] = find(label_refined == m);
                len2 = max(row2) - min(row2) + 1;
                breadth2 = max(col2) - min(col2) + 1;
                target6  = uint8(zeros([len2 breadth2]));
                target7  = uint8(zeros([len2 breadth2]));
                target8  = zeros([len2 breadth2]);
                target9  = zeros([len2 breadth2]);
                target10 = zeros([len2 breadth2]);
                sy2 = min(col2) - 1;
                sx2 = min(row2) - 1;

                for n = 1:length(row2)
                    x2 = row2(n) - sx2;
                    y2 = col2(n) - sy2;
                    target6(x2, y2)  = target1(row2(n), col2(n)); % Original Bright-field image
                    target7(x2, y2)  = target2(row2(n), col2(n)); % Original Green
                    target8(x2, y2)  = target3(row2(n), col2(n)); % Bright-field Mask
                    target9(x2, y2)  = target4(row2(n), col2(n)); % Nucleus Mask
                    target10(x2, y2) = target5(row2(n), col2(n)); % Green Particle Mask
                end

                label_Nucleus = bwlabel(target9);
                Nucleus_count = max(label_Nucleus(:));
                label_GreenMask_multi = bwlabel(target10);
                num_particles_multi = max(label_GreenMask_multi(:));

                if Nucleus_count == 1
                    nuc_erode = imerode(target9,strel('diamond',3));

                    if num_particles_multi == 1 
                        imwrite(target6, fullfile(CellBFSingle_path, sprintf('%d.tif', ImgCellSingle_count + 1)));
                        imwrite(target7, fullfile(CellGreenSingle_path, sprintf('%d.tif', ImgCellSingle_count + 1))); % Optional
                        ImgCellSingle_count = ImgCellSingle_count + 1;
                    elseif num_particles_multi > 1
                        imwrite(target6, fullfile(CellBFMulti_path, sprintf('%d.tif', ImgCellMulti_count + 1)));
                        imwrite(target7, fullfile(CellGreenMulti_path, sprintf('%d.tif', ImgCellMulti_count + 1))); % Optional
                        ImgCellMulti_count = ImgCellMulti_count + 1;
                    else
                        imwrite(target6, fullfile(CellBFNone_path, sprintf('%d.tif', ImgCellNone_count + 1)));
                        imwrite(target7, fullfile(CellGreenNone_path, sprintf('%d.tif', ImgCellNone_count + 1))); % Optional
                        ImgCellNone_count = ImgCellNone_count + 1;
                    end
                end
            end
        end
    end
end

%% =========================================================================
% 
%                   ResNet Classification Model
%
% ==========================================================================

% =================================================================
%  PART 1: TRAIN THE 3-CATEGORY RESNET MODEL
% =================================================================

%% === 1. Load and Prepare the 3-Category Data ===

trainingData_path = 'D:/Image Process DNA/Training (ResNet)/PNI + DNA';

% --- Category 1: Single Particle ---
single_folder = fullfile(trainingData_path, 'Single Particle', 'BF');
single_files = dir(fullfile(single_folder, '*.tif'));
single_paths = fullfile(single_folder, {single_files.name}');
single_labels = repelem(categorical({'SingleParticle'}), numel(single_paths), 1);

% --- Category 2: Multiple Particle ---
multi_folder = fullfile(trainingData_path, 'Multiple Particle', 'BF');
multi_files = dir(fullfile(multi_folder, '*.tif'));
multi_paths = fullfile(multi_folder, {multi_files.name}');
multi_labels = repelem(categorical({'MultiParticle'}), numel(multi_paths), 1);

% --- Category 3: No Particle ---
none_folder = fullfile(trainingData_path, 'No Particle', 'BF');
none_files = dir(fullfile(none_folder, '*.tif'));
none_paths = fullfile(none_folder, {none_files.name}');
none_labels = repelem(categorical({'NoParticle'}), numel(none_paths), 1);

% --- Combine all data ---
allFilePaths = [single_paths; multi_paths; none_paths];
allLabels = [single_labels; multi_labels; none_labels];

% --- Create the imageDatastore ---
imds = imageDatastore(allFilePaths, 'Labels', allLabels);
imds = shuffle(imds);
[imdsTrain, imdsValidation] = splitEachLabel(imds, 0.8, 'randomized');

%% === 2. Load and Adapt the Pre-Trained ResNet-18 Model ===

net = resnet18;
lgraph = layerGraph(net);
inputSize = net.Layers(1).InputSize;

% --- Adapt the first layer for grayscale input ---
layers = lgraph.Layers;
firstConvLayer = layers(2); 
inputWeights = firstConvLayer.Weights;
avgWeights = mean(inputWeights, 3);
newWeights = cat(3, avgWeights, avgWeights, avgWeights);
newFirstConvLayer = convolution2dLayer(firstConvLayer.FilterSize, ...
    firstConvLayer.NumFilters, 'Name', firstConvLayer.Name, ...
    'Weights', newWeights, 'Bias', firstConvLayer.Bias, ...
    'Stride', firstConvLayer.Stride, 'Padding', firstConvLayer.Padding);
lgraph = replaceLayer(lgraph, firstConvLayer.Name, newFirstConvLayer);


%% === 3. Replace the Final Layers for Our 3-Class Task ===

numClasses = numel(categories(imdsTrain.Labels));
classNames = categories(imdsTrain.Labels);
pref = containers.Map( ...
    {'SingleParticle','MultiParticle','NoParticle'}, ...
    [4.5,              3.5,             1.0]); 

classWeights = zeros(numel(classNames),1);
for k = 1:numel(classNames)
    classWeights(k) = pref(char(classNames{k}));
end

learnableLayerName   = 'fc1000';
classOutputLayerName = 'ClassificationLayer_predictions';

newLearnableLayer = fullyConnectedLayer(numClasses, ...
    'Name','new_fc', 'WeightLearnRateFactor',10, 'BiasLearnRateFactor',10);

newClassLayer = classificationLayer( ...
    'Name','new_classoutput', ...
    'Classes',classNames, ...
    'ClassWeights',classWeights);

lgraph = replaceLayer(lgraph, learnableLayerName, newLearnableLayer);
lgraph = replaceLayer(lgraph, classOutputLayerName, newClassLayer);


%% === 4. Set up Training Options and Train the Network ===

imageAugmenter = imageDataAugmenter( ...
    'RandXReflection', true, ...
    'RandRotation',   [-10, 10], ...
    'RandScale',      [0.92, 1.08], ...
    'RandXTranslation',[-8, 8], ...
    'RandYTranslation',[-8, 8]);
augimdsTrain = augmentedImageDatastore(inputSize(1:2), imdsTrain, ...
    'DataAugmentation', imageAugmenter, 'ColorPreprocessing','gray2rgb');
augimdsValidation = augmentedImageDatastore(inputSize(1:2), imdsValidation, ...
    'ColorPreprocessing','gray2rgb');

valFreq = max(5, floor(numel(imdsTrain.Files)/ (10*32)));

options = trainingOptions('adam', ...
    'MaxEpochs', 50, ...
    'InitialLearnRate', 3e-4, ... 
    'MiniBatchSize', 32, ...
    'Shuffle','every-epoch', ...
    'ValidationData', augimdsValidation, ...
    'ValidationFrequency', valFreq, ...
    'ValidationPatience', 12, ...
    'OutputNetwork','best-validation-loss', ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod', 15, ...
    'LearnRateDropFactor', 0.2, ...
    'GradientThresholdMethod','l2norm', 'GradientThreshold',5, ...
    'L2Regularization', 5e-5, ...  
    'Plots','training-progress');

% --- TRAIN THE NETWORK ---
[myTrainedNet, trainInfo] = trainNetwork(augimdsTrain, lgraph, options);

% === Convert iteration-level logs into epoch-level averages ===
numEpochs = options.MaxEpochs;
numIterations = numel(trainInfo.TrainingLoss);
iterationsPerEpoch = floor(numIterations / numEpochs);

trainLoss_epoch = zeros(1,numEpochs);
trainAcc_epoch  = zeros(1,numEpochs);

for e = 1:numEpochs
    idxStart = (e-1)*iterationsPerEpoch + 1;
    if e < numEpochs
        idxEnd = e*iterationsPerEpoch;
    else
        idxEnd = numIterations;
    end
    trainLoss_epoch(e) = mean(trainInfo.TrainingLoss(idxStart:idxEnd));
    trainAcc_epoch(e)  = mean(trainInfo.TrainingAccuracy(idxStart:idxEnd));
end

% --- Epoch indices for plotting ---
epochs = 1:numEpochs;
valEpochs = linspace(1,numEpochs,numel(trainInfo.ValidationLoss)); % spread validation points

figure;
plot(epochs, trainLoss_epoch, '-k','LineWidth',1.5); hold on;
plot(valEpochs, trainInfo.ValidationLoss, '-r','LineWidth',1.5);
xlabel('Epoch'); ylabel('Loss');
legend('Training','Validation','Location','best'); legend boxoff;
title('Training vs Validation Loss');
grid on;

figure;
plot(epochs, trainAcc_epoch, '-k','LineWidth',1.5); hold on;
plot(valEpochs, trainInfo.ValidationAccuracy, '-r','LineWidth',1.5);
xlabel('Epoch'); ylabel('Accuracy (%)');
legend('Training','Validation','Location','best'); legend boxoff;
title('Training vs Validation Accuracy');
grid on;

%% === 5. Evaluate Model Performance ===

Y_true = imdsValidation.Labels;
Y_pred = classify(myTrainedNet, augimdsValidation);

figure;
cm_chart = confusionchart(Y_true, Y_pred);
cm_chart.Title = 'Confusion Matrix for 3-Class BF-Only Model';
cm_chart.RowSummary = 'row-normalized'; 
cm_chart.ColumnSummary = 'column-normalized';

% accuracy = sum(Y_pred == Y_true) / numel(Y_true);
% fprintf('Overall Validation Accuracy: %.2f%%\n', accuracy * 100);
% 
% save('ResNetCellClassifier_3Class_Final.mat', 'myTrainedNet', 'trainInfo');
% 
% %% =================================================================
% %  PART 2: APPLY THE TRAINED MODEL
% % =================================================================
% 
% %% === 6. Load the Trained Model ===
% 
% load('ResNetCellClassifier_3Class_Final.mat');
% inputSize = myTrainedNet.Layers(1).InputSize;
% 
% %% === 7. Point to the FOLDER of New Images to Classify ===
% 
% Apply_path = 'D:/Image Process DNA/Training (imfindcircle)/PNI + DNA/BF';
% Apply_files = get_files_numerically(Apply_path);
% 
% %% === 8. Classify Each Image and Store the Result ===
% 
% imdsApply = imageDatastore(Apply_path);
% augimdsApply = augmentedImageDatastore(inputSize(1:2), imdsApply, 'ColorPreprocessing', 'gray2rgb');
% [predictedLabels, probabilities] = classify(myTrainedNet, augimdsApply);
% resultsTable = table();
% resultsTable.FileName = {Apply_files.name}';
% resultsTable.PredictedLabel = predictedLabels;
% resultsTable.Confidence = max(probabilities, [], 2);
% 
% %% === 9. Final Summary Report ===
% 
% fprintf('\n--- Final Summary Report ---\n');
% disp('--- First 10 Classification Results ---');
% disp(head(resultsTable, 10));
% 
% total_cells_analyzed = height(resultsTable);
% count_single = sum(resultsTable.PredictedLabel == 'SingleParticle');
% count_multi = sum(resultsTable.PredictedLabel == 'MultiParticle');
% count_none = sum(resultsTable.PredictedLabel == 'NoParticle');
% 
% fprintf('\n=========================================\n');
% fprintf('         FINAL DETECTION SUMMARY\n');
% fprintf('=========================================\n');
% fprintf('Total cells analyzed: %d\n', total_cells_analyzed);
% fprintf('-----------------------------------------\n');
% fprintf('Cells classified as "SingleParticle": %d\n', count_single);
% fprintf('Cells classified as "MultiParticle":  %d\n', count_multi);
% fprintf('Cells classified as "NoParticle":     %d\n', count_none);
% fprintf('=========================================\n');
% 

T = table( ...
    (1:numel(trainLoss_epoch)).', ...
    trainLoss_epoch(:), ...
    trainAcc_epoch(:), ...
    linspace(1,numel(trainLoss_epoch),numel(trainInfo.ValidationLoss)).', ...
    trainInfo.ValidationLoss(:), ...
    trainInfo.ValidationAccuracy(:), ...
    'VariableNames', {'Epoch','TrainLoss','TrainAccuracy','ValEpoch','ValLoss','ValAccuracy'});

writetable(T, 'training_log.xlsx');

% data = readtable('training_log.xlsx');
% figure;
% plot(data.Epoch, data.ValidationAccuracy, 'r-', 'LineWidth', 2);
% hold on;
% plot(data.Epoch, data.ValidationLoss, 'b--', 'LineWidth', 2);
% xlabel('Epoch');
% ylabel('Value');
% legend('Validation Accuracy', 'Validation Loss');
% title('Validation Performance');