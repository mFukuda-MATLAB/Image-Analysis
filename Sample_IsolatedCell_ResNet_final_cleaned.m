% =========================================================================
%
% ResNet-18 Application Pipeline: Single-Particle Selection & FITC Measurement
%
% Purpose:
%   Apply a trained 3-class ResNet-18 model to segmented brightfield ROIs
%   to automatically select cells containing exactly one phagocytosed
%   microparticle, then confirm with circle detection and quantify FITC
%   intensity within the particle region.
%
% Overview of Steps:
%   1) Load phase (BF), DAPI, and FITC images (per-sample folders).
%   2) Segment single-cell ROIs using DAPI-guided masks (isolates single
%      nuclei and refines cell boundaries in BF).
%   3) Classify each ROI with a pretrained ResNet-18 (Single/Multi/None);
%      keep only predicted SingleParticle images (optional confidence gate).
%   4) Confirm a single particle by matched circle detection:
%        - BF “halo” (bright ring) vs. LoG-based “core” (dark disk)
%        - geometric proximity match → final (center, radius)
%   5) Measure FITC intensity within the detected particle mask
%      (optionally dilated by a fraction of the radius) and save results.
%
% Inputs (8-bit TIFFs; pre-split by ImageJ/Fiji):
%   - Brightfield (Phase):  ParticleType_##.tif
%   - DAPI (Blue):         ParticleType_##.tif
%   - FITC (Green):        ParticleType_##.tif
%   Directory structure is discovered as Particle_Type / Sample / {Phase,Blue,Green}.
%
% Outputs:
%   - Single-cell BF crops (and matched Green/NoNuc versions) for each sample
%   - ResNet-positive images copied to:  ResNet Results/{BF, Green, No Nucleus}
%   - Final confirmed single-particle images copied to:
%       • Final_Single_Particles (BF)
%       • Final_Single_Particles (Green)
%   - Table of per-image FITC mean intensity within the particle mask
%     (one row per confirmed single-particle cell)
%
% Model Notes:
%   - Backbone: pretrained ResNet-18 (ImageNet)
%   - Input handling: grayscale BF resized to 224×224 and expanded to 3-ch
%   - Positive class = 'SingleParticle' (optional confidence threshold)
%
% Key Parameters (this script):
%   - imfindcircles ranges: BF halos [4 12], LoG cores [4 12]
%   - Minimum particle area, intensity/solidity thresholds
%   - Dilation fraction for measurement mask (e.g., 0.20)
%
% Requirements:
%   - MATLAB R2023a or later
%   - Image Processing Toolbox
%   - Deep Learning Toolbox (+ pretrained ResNet-18 support package)
%   - ImageJ/Fiji (channel decomposition)
%
% Reproducibility (recommended):
%   - rng(42,'twister'); record GPU/CUDA info and MATLAB version
% =========================================================================

clc; clear; warning('off','all'); close all

% Directory (Example directory)
main_dir = 'D:/Image Process DNA/DNA Particle';

Particle_Type_folder_info = dir(main_dir);
Particle_Type = {Particle_Type_folder_info([Particle_Type_folder_info.isdir] & ...
    ~ismember({Particle_Type_folder_info.name}, {'.', '..'})).name};

for a = 1:length(Particle_Type)

    Current_folder = fullfile(main_dir, Particle_Type{a});
    Sample_folder_info = dir(Current_folder);
    Samples = {Sample_folder_info([Sample_folder_info.isdir] & ...
    ~ismember({Sample_folder_info.name}, {'.', '..'})).name};

    master_intensity_results = table();

    for b = 1:length(Samples)
        % Directory (Example directory & Pathes)
        SingleCell_path   = fullfile(main_dir, Particle_Type{a}, Samples{b} ,'Single Cells (BF)');         mkdir(SingleCell_path);
        Green_path        = fullfile(main_dir, Particle_Type{a}, Samples{b} ,'Single Cells (Green)');      mkdir(Green_path);
        NoNuc_path        = fullfile(main_dir, Particle_Type{a}, Samples{b} ,'Single Cells (No Nucleus)'); mkdir(NoNuc_path);

        ResNetResult_path   = fullfile(main_dir, Particle_Type{a}, Samples{b} ,'ResNet Results');
        BF_pos_path         = fullfile(ResNetResult_path, 'BF');         mkdir(BF_pos_path);
        Green_pos_path      = fullfile(ResNetResult_path, 'Green');      mkdir(Green_pos_path);
        NoNuc_pos_path      = fullfile(ResNetResult_path, 'No Nucleus'); mkdir(NoNuc_pos_path);
    
        ImageCount = 0;
    
        phase_dir  = fullfile(main_dir, Particle_Type{a}, Samples{b}, 'Phase');
        blue_dir   = fullfile(main_dir, Particle_Type{a}, Samples{b}, 'Blue');
        green_dir  = fullfile(main_dir, Particle_Type{a}, Samples{b}, 'Green');

        phase_files = dir(fullfile(phase_dir, '*.tif'));
        nums = zeros(1, numel(phase_files));

        for aa = 1:numel(phase_files)
            name = phase_files(aa).name;
            match = regexp(name, '\d+', 'match');
            if ~isempty(match)
                nums(aa) = str2double(match{end});  % Use last number if multiple found
            else
                nums(aa) = NaN;
            end
        end

% ==========================================================================
%
%                       PART 1: Cell Segmentation
%
% ==========================================================================

        % Sort by numeric value and apply to file list
        [~, idx] = sort(nums);
        phase_files = phase_files(idx);

        contrast_list = [];
        noise_list = [];

        for c = 1:length(phase_files)
            img = im2double(imread(fullfile(phase_dir, phase_files(c).name)));
            contrast_list(end+1) = std2(img);
            noise_list(end+1) = mean2(stdfilt(img));
        end

        % Threshold for Processing Images
        Contrast_threshold = prctile(contrast_list, 20);
        Noise_threshold = prctile(noise_list, 75);

        Img_count_ParCell = 0;

        for d = 1:length(phase_files)

            phase_filename = phase_files(d).name;
            blue_filename  = strrep(phase_filename, 'Phase', 'Blue');
            green_filename  = strrep(phase_filename, 'Phase', 'Green');
            Img_1 = fullfile(phase_dir, phase_filename);
            Img_2 = fullfile(blue_dir, blue_filename);
            Img_3 = fullfile(green_dir, green_filename);

            if ~isfile(Img_1) || ~isfile(Img_2) || ~isfile(Img_3)
                continue;
            end

        %% Reading images

            BF = imread(Img_1);
            Blue = imread(Img_2);
            Green = imread(Img_3);

            [H, W] = size(BF);
            clean_mask = false(H, W);

        %% Parameters (Convolution)

             % For BF
            CF1 = fspecial('sobel'); Trans_CF1 = CF1';
            Conv_Threshold_BF = 0.2;

            % For Nucleus
            CF2 = fspecial('gaussian'); Trans_CF2 = CF2';
            Conv_Threshold_Nuc = 1.1;

            % Threshold for Processing Images
            BW_noise_threshold = 0.15;

        %% Image Processing (Convolution Mask)

            % Processing Bright-field Images
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

            % Making Bright-field mask
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

            % figure()
            % imshow(imoverlay(BF,bwperim(mask_clean)))

            % Processing Blue Images
            I_Nuc = im2double(Blue);
            I_med_Nuc = imgaussfilt(I_Nuc,0.8);
            edges_Nuc = sqrt((conv2(I_med_Nuc, CF2, 'same')).^2 + (conv2(I_med_Nuc, Trans_CF2, 'same')).^2);
            T_Nuc = graythresh(edges_Nuc);
            bw_Nuc = imbinarize(edges_Nuc, T_Nuc * Conv_Threshold_Nuc);
            bw_Nuc_clean = imclose(bw_Nuc,strel('disk',10));
            bw_Nuc_clean = imdilate(bw_Nuc_clean, strel('disk', 8));
            bw_Nuc_clean = imopen(bw_Nuc_clean, ones(5,5));
            bw_Nuc_clean = imfill(bw_Nuc_clean,'holes');
            bw_Nuc_clean = imerode(bw_Nuc_clean,strel('diamond',10));

            % Processing Green Images
            I_green = imbinarize(Green);
            D = bwdist(~I_green);

            bw_Green_clean = imextendedmax(D,2);
            bw_Green_clean = imdilate(bw_Green_clean, strel('disk',10));
            bw_Green_clean = imclearborder(bw_Green_clean);
            bw_Par = imerode(bw_Green_clean, strel('disk', 5));

            label_green = bwlabel(bw_Par);
            label_Mask = bwlabel(mask_clean);
            label_Nuc = bwlabel(bw_Nuc_clean);
            n_cells = max(label_Mask(:));

            % figure()
            % imshow(imoverlay(BF,bwperim(mask_clean)));

            %% Checking the efficiency of cell segmentation

            % label_mask = bwlabel(bw_BF_clean);
            % cent_num = regionprops(label_mask,'Centroid');
            % 
            % for bb = 1:numel(cent_num)
            %     text(cent_num(bb).Centroid(1),cent_num(bb).Centroid(2),num2str(bb),'color','r')
            % end
            % 
            % figure()
            % imshow(imoverlay(BF,bwperim(bw_BF_clean)))
            % figure();
            % imshow(imoverlay(imoverlay(imoverlay(mat2gray(BF),mask_clean,'blue'),bw_Nuc_clean,'red'),bw_Par,'green'));
            % pause(1); close;

         %% Cell Segmentation

            for j = 1:n_cells
                cell_mask = (label_Mask == j);
                [row1, col1] = find(cell_mask);
                len = max(row1) - min(row1) + 1;
                breadth = max(col1) - min(col1) + 1;
                target1 = uint8(zeros([len breadth]));
                target2 = zeros([len breadth]);
                target3 = zeros([len breadth]);
                sy = min(col1) - 1;
                sx = min(row1) - 1;

                    for k = 1:length(row1)
                        x = row1(k) - sx;
                        y = col1(k) - sy;
                        target1(x, y) = BF(row1(k), col1(k));
                        target2(x, y) = bw_Nuc_clean(row1(k), col1(k));
                        target3(x, y) = mask_clean(row1(k), col1(k));
                    end

                 label_Nuc_single = bwlabel(target2);

            % Isolated Cells
                if max(label_Nuc_single(:)) == 1

                    clean_mask(sx + (1:len), sy + (1:breadth)) = clean_mask(sx + (1:len), sy + (1:breadth)) | (target3 > 0);

            % Aggregated Cells
                elseif max(label_Nuc_single(:)) > 1

                    I = im2double(target1);
                    I = adapthisteq(I, 'ClipLimit', 0.05, 'Distribution', 'exponential', ...
                                    'Alpha', 10, 'NBins', 200);
                    I = imsharpen(I, 'Radius', 1.5, 'Amount', 2.5, 'Threshold', 0.01);
                    I = imgaussfilt(I, 1.8);

                    G = imgradient(I, 'sobel');
                    G = imbilatfilt(G,std2(G)^2, 5);          
                    G = imadjust(G, stretchlim(G), []);
                    G = imcomplement(G);           

                    nuclei     = imerode(target2 > 0, strel('disk', 3)); 
                    mask_crop  = imerode(target3 > 0, strel('disk', 3)); 

                    G2 = imimposemin(G, nuclei);
                    L  = watershed(G2);
                    L(~mask_crop) = 0;
                    L  = bwlabel(L > 0); 

                    for m = 1:max(L(:))
                        region = (L == m);

                        [row2, col2] = find(region);
                        minR = min(row2); maxR = max(row2);
                        minC = min(col2); maxC = max(col2);

                        mask_crop = region(minR:maxR, minC:maxC);
                        len2 = maxR - minR + 1;
                        breadth2 = maxC - minC + 1;

                        sx2 = minR - 1;
                        sy2 = minC - 1;

                        x_start = sx + sx2 + 1;
                        x_end   = x_start + len2 - 1;
                        y_start = sy + sy2 + 1;
                        y_end   = y_start + breadth2 - 1;

                        % Safety check to avoid exceeding image bounds
                        if x_end > H || y_end > W || x_start < 1 || y_start < 1
                            continue;
                        end

                        % Place into full-size mask
                        clean_mask(x_start:x_end, y_start:y_end) = ...
                        clean_mask(x_start:x_end, y_start:y_end) | mask_crop;
                    end
                end          
            end

            [H, W] = size(clean_mask);
            labeled = bwlabel(clean_mask);
            props = regionprops(labeled, 'PixelIdxList');
            clean_mask_filtered = false(H, W);

            % === Set range for boundary % to remove cells on border===
            min_frac = 0.10;   % 10% or more
            max_frac = 0.70;   % but not more than 70%

            for i = 1:length(props)
                obj_mask = false(H, W);
                obj_mask(props(i).PixelIdxList) = true;

                % Get boundary
                B = bwboundaries(obj_mask);
                if isempty(B), continue; end
                boundary = B{1};  % [row, col]

                num_total = size(boundary, 1);
                on_border = boundary(:,1) == 1 | boundary(:,1) == H | ...
                            boundary(:,2) == 1 | boundary(:,2) == W;
                num_border = sum(on_border);

                frac = num_border / num_total;

                % Remove if within "truncation suspicious" range
                if frac < min_frac || frac > max_frac
                    clean_mask_filtered = clean_mask_filtered | obj_mask;
                end
            end

            clean_mask = clean_mask_filtered;

            label_clean = bwlabel(clean_mask);
            refined_mask = false(size(clean_mask));

            for l = 1:max(label_clean(:))
                region = (label_clean == l);
                if nnz(region) < 100
                    continue;
                end
                initMask = imerode(region, strel('disk', 3));
                refined = activecontour(I_med_BF, initMask, 10, 'edge');
                refined = imdilate(refined, strel('disk', 1));
                refined_mask = refined_mask | refined;
            end

            clean_mask = refined_mask;
            % figure()
            % imshow(imoverlay(BF, bwperim(clean_mask)))

            label_final = bwlabel(clean_mask);

            for o = 1:max(label_final(:))
                [row3, col3] = find(label_final == o);
                len3 = max(row3) - min(row3) + 1;
                breadth3 = max(col3) - min(col3) + 1;
                target5 = uint8(zeros([len3 breadth3]));
                target6 = zeros([len3 breadth3]);
                target7 = uint8(zeros([len3 breadth3]));
                sy3 = min(col3) - 1;
                sx3 = min(row3) - 1;

                for p = 1:length(row3)
                    x3 = row3(p) - sx3;
                    y3 = col3(p) - sy3;
                    target5(x3, y3) = BF(row3(p), col3(p)); % Original Bright-field image
                    target6(x3, y3) = bw_Nuc_clean(row3(p), col3(p)); % Nuclear Mask
                    target7(x3, y3) = Green(row3(p), col3(p)); % Green
                end

                if max(target6(:)) == 1
                imwrite(target5, fullfile(SingleCell_path, sprintf('%d.tif', ImageCount + 1)));
                imwrite(target7, fullfile(Green_path, sprintf('%d.tif', ImageCount + 1)));
                imwrite(target5 .* uint8(~target6), fullfile(NoNuc_path, sprintf('%d.tif', ImageCount + 1)));
                ImageCount = ImageCount + 1;
                end
            end
        end

%% =================================================================
% PART 2: APPLY THE TRAINED 3-CATEGORY MODEL
% ====================================================================

        % === 1. Load the Trained Model and Get its Input Size ===

        % Load the .mat file created by your training script
        load('ResNetCellClassifier_3Class_Re.mat');
        inputSize = myTrainedNet.Layers(1).InputSize;

        % === 2. Point to the FOLDER of New Images to Classify ===

        Apply_path = fullfile(SingleCell_path);
        Apply_files = get_files_numerically(Apply_path);

        % === 3. Classify Each Image and Store the Result ===

        % Pre-allocate a table to hold the results for clarity
        numImages = length(Apply_files);
        resultsTable = table('Size', [numImages, 3], ...
                             'VariableTypes', {'string', 'categorical', 'double'}, ...
                             'VariableNames', {'FileName', 'PredictedLabel', 'Confidence'});

        positiveClassName = 'SingleParticle';
        confidence_threshold = 0.75; % You can tune this value

        for i = 1:numImages
            % --- Load the next cell image ---
            imgPath = fullfile(Apply_path, Apply_files(i).name);
            newCellImage = imread(imgPath);

            green_imgPath = fullfile(Green_path, Apply_files(i).name);
            NoNuc_imgPath = fullfile(NoNuc_path, Apply_files(i).name);

            % 1. Resize the image
            imgResized = imresize(newCellImage, inputSize(1:2));
            % 2. If the image is grayscale, convert it to 3-channel RGB
            if size(imgResized, 3) == 1
                imgResized = cat(3, imgResized, imgResized, imgResized);
            end

            % --- Classify the image ---
            [predictedLabel, probability] = classify(myTrainedNet, imgResized);
            predictedLabel = string(predictedLabel);

            % --- Filter Out Low-Confidence Positives ---
            if strcmp(predictedLabel, positiveClassName)
               if max(probability) < confidence_threshold
                    predictedLabel = 'NoParticle'; 
               end
            end

            % --- Store the results in our table ---
            resultsTable.PredictedLabel(i) = predictedLabel;
            resultsTable.FileName(i) = Apply_files(i).name;
            resultsTable.Confidence(i) = max(probability);

            if strcmp(predictedLabel, positiveClassName)
                copyfile(imgPath, BF_pos_path);
                copyfile(green_imgPath, Green_pos_path);
                copyfile(NoNuc_imgPath, NoNuc_pos_path);
            end
        end

        % === 4. Final Summary Report ===
        excel_save_path = fullfile(ResNetResult_path, 'resultsTable.xlsx');
        writetable(resultsTable, excel_save_path);

        % --- Calculate the final summary statistics ---
        total_cells_analyzed = height(resultsTable);
        count_single = sum(strcmp(cellstr(resultsTable.PredictedLabel), 'SingleParticle'));
        count_multi = sum(strcmp(cellstr(resultsTable.PredictedLabel), 'MultiParticle'));
        count_none = sum(strcmp(cellstr(resultsTable.PredictedLabel), 'NoParticle'));

        fprintf('\n=========================================\n');
        fprintf('         FINAL DETECTION SUMMARY\n');
        fprintf('=========================================\n');
        fprintf('Total cells analyzed: %d\n', total_cells_analyzed);
        fprintf('-----------------------------------------\n');
        fprintf('Cells classified as "SingleParticle": %d\n', count_single);
        fprintf('Cells classified as "MultiParticle":  %d\n', count_multi);
        fprintf('Cells classified as "NoParticle":     %d\n', count_none);
        fprintf('=========================================\n');

        if total_cells_analyzed > 0
            percentage_with_particles = (count_single / total_cells_analyzed) * 100;
            fprintf('Percentage of cells with single particles: %.2f%%\n', percentage_with_particles);
        end
        fprintf('=========================================\n');

%% =================================================================
%  PART 3: FINAL FILTERING AND INTENSITY MEASUREMENT (MATCHED METHOD)
% ==================================================================

        BF_NoNuc_files = get_files_numerically(NoNuc_pos_path);

        SINGLE_PARTICLE_BF_PATH = fullfile(main_dir, Particle_Type{a}, Samples{b}, 'Final_Single_Particles (BF)');
        mkdir(SINGLE_PARTICLE_BF_PATH);
        SINGLE_PARTICLE_Green_PATH = fullfile(main_dir, Particle_Type{a}, Samples{b}, 'Final_Single_Particles (Green)');
        mkdir(SINGLE_PARTICLE_Green_PATH);

        % --- Parameters for Filtering ---
        PARTICLE_SIGMA = 3;
        FILTER_SIZE = 2 * ceil(PARTICLE_SIGMA*3) + 1;
        MIN_PARTICLE_AREA = 130;

        % --- Parameters for Measurement ---
        DILATION_FRACTION = 0.20;   % Same as the red-mask code

        sample_intensity_results = table();
        num_saved = 0;

        for z = 1:length(BF_NoNuc_files)

            current_filename = BF_NoNuc_files(z).name;

            bf_image_path       = fullfile(BF_pos_path, BF_NoNuc_files(z).name);
            green_image_path    = fullfile(Green_pos_path, BF_NoNuc_files(z).name);
            bf_NoNuc_image_path = fullfile(NoNuc_pos_path, BF_NoNuc_files(z).name);

            bf_image_data = imread(bf_image_path);
            Green_image_data = imread(green_image_path);
            bf_NoNuc_image_data = imread(bf_NoNuc_image_path);

            % --- Detect cores and halos (same as before) ---
            bf_double = mat2gray(bf_NoNuc_image_data);
            log_filter = fspecial('gaussian', FILTER_SIZE, PARTICLE_SIGMA);
            log_filtered_image = imfilter(bf_double, log_filter, 'replicate');
            log_response = mat2gray(log_filtered_image);
            processed_image = imcomplement(log_response);
            bright_processed_image = imadjust(processed_image, [0.1 0.6], [0.0 1.0]);

            [centers_log, radii_log, ~] = imfindcircles(bright_processed_image, [4 12], ...
                'Sensitivity', 0.85,'EdgeThreshold', 0.3, 'ObjectPolarity', 'dark');

            bf_brightened = imadjust(bf_double, stretchlim(bf_double, [0.10 0.90]));
            bf_adjusted = adapthisteq(bf_brightened, 'ClipLimit', 0.05, 'Distribution','rayleigh');
            bf_smoothed = imgaussfilt(bf_adjusted, 1.5);

            [centers_bf, radii_bf, ~] = imfindcircles(bf_smoothed, [4 12], ...
                'Sensitivity', 0.83,'EdgeThreshold', 0.01, 'ObjectPolarity', 'bright');

            % --- Select valid halos and cores ---
            MIN_HALO_INTENSITY = 0.6;
            MIN_HALO_SOLIDITY = 0.8;
            valid_halos_centers = [];
            valid_halos_radii = [];
            for w = 1:size(centers_bf, 1)
                [X,Y] = meshgrid(1:size(bf_smoothed,2), 1:size(bf_smoothed,1));
                circle_mask = (X - centers_bf(w,1)).^2 + (Y - centers_bf(w,2)).^2 <= radii_bf(w)^2;

                props = regionprops(circle_mask, bf_smoothed, 'Area', 'MeanIntensity', 'Solidity');
                if isempty(props), continue; end

                if (props.MeanIntensity >= MIN_HALO_INTENSITY) && (props.Solidity >= MIN_HALO_SOLIDITY) && ...
                        (props.Area >= MIN_PARTICLE_AREA)
                    valid_halos_centers = [valid_halos_centers; centers_bf(w,:)];
                    valid_halos_radii = [valid_halos_radii; radii_bf(w)];
                end
            end

            MAX_CORE_INTENSITY = 0.4;
            MIN_CORE_SOLIDITY = 0.90;
            valid_cores_centers = [];
            valid_cores_radii = [];
            for w = 1:size(centers_log, 1)
                [X,Y] = meshgrid(1:size(bright_processed_image,2), 1:size(bright_processed_image,1));
                circle_mask = (X - centers_log(w,1)).^2 + (Y - centers_log(w,2)).^2 <= radii_log(w)^2;

                props = regionprops(circle_mask, bright_processed_image, 'Area', 'MeanIntensity', 'Solidity');
                if isempty(props), continue; end

                if (props.MeanIntensity <= MAX_CORE_INTENSITY) && (props.Solidity >= MIN_CORE_SOLIDITY) && ...
                        (props.Area >= MIN_PARTICLE_AREA)
                    valid_cores_centers = [valid_cores_centers; centers_log(w,:)];
                    valid_cores_radii = [valid_cores_radii; radii_log(w)];
                end
            end

            % --- Match halos and cores ---
            final_centers = [];
            final_radii = [];
            if ~isempty(valid_halos_centers) && ~isempty(valid_cores_centers)
                match_threshold = 10.0;
                for i = 1:size(valid_cores_centers, 1)
                    distances = pdist2(valid_cores_centers(i,:), valid_halos_centers);
                    if min(distances) < match_threshold
                        final_centers = [final_centers; valid_cores_centers(i,:)];
                        final_radii = [final_radii; valid_cores_radii(i)];
                    end
                end
            end

            % --- Final selection and measurement ---
            if size(final_centers, 1) == 1

                saveBF_path = fullfile(SINGLE_PARTICLE_BF_PATH, current_filename);
                saveGreen_path = fullfile(SINGLE_PARTICLE_Green_PATH, current_filename);
                % imwrite(bf_image_data, saveBF_path);
                % imwrite(Green_image_data, saveGreen_path);
                num_saved = num_saved + 1;

                % === Match the same logic as the red-mask method ===
                [N, M] = meshgrid(1:size(Green_image_data,2), 1:size(Green_image_data,1));
                redMask = (N - final_centers(1,1)).^2 + (M - final_centers(1,2)).^2 <= final_radii^2;

                % Dilate slightly by 20% of radius
                dilation_radius = round(DILATION_FRACTION * (final_radii / 2));
                se = strel('disk', dilation_radius);
                padded_redMask = imdilate(redMask, se);

figure;

% Left: original brightfield image
subplot(1,2,1);
imshow(bf_image_data, []);

% Right: overlay with final detection
subplot(1,2,2);
imshow(bf_image_data, []); hold on;
viscircles(final_centers(1,:), double(final_radii), 'Color', 'r', 'LineWidth', 2);
hold off;
                

                % === Calculate mean green intensity ===
                green_intensity_img_raw = im2double(Green_image_data);
                green_intensities = green_intensity_img_raw(padded_redMask);

                if isempty(green_intensities)
                    mean_particle_intensity = NaN;
                else
                    mean_particle_intensity = mean(green_intensities);
                end

                % === Store results ===
                new_row = table();
                new_row.SampleName = Samples{b};
                new_row.FileName = string(BF_NoNuc_files(z).name);
                new_row.RawParticleIntensity = mean_particle_intensity;

                sample_intensity_results = [sample_intensity_results; new_row];
            end
        end

        % === Combine and export ===
        % master_intensity_results = [master_intensity_results; sample_intensity_results];
        % excel_filename = sprintf('%s_intensity_measurements_rev02.xlsx', Particle_Type{a});
        % excel_save_path = fullfile(main_dir, Particle_Type{a}, excel_filename);
        % writetable(master_intensity_results, excel_save_path);
    end
end


