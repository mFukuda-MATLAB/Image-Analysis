%% =========================================================================
%
% His/DNA +LysoView Segmentation and Measurement (with PAH-AF594)
%
% Mask images:
%       Binary segmentation masks generated from brightfield (cell boundary),
%       DAPI (nucleus), and FITC (particle) channels.
%       These masks are used to identify ROIs and classify cells into
%       single-particle, multi-particle, or no-particle categories.
%
% File inputs: (decomposed TIFF images into separate channels using ImageJ)
%       FITC images  (8 bit tifs with format ParticleType_##.tif)
%       DAPI images  (8 bit tifs with format ParticleType_##.tif)
%       Phase (Bright Field) images (8 bit tifs with format ParticleType_##.tif)
%
% ResNet Training: (3-classification model)
%       The model was trained on cell images containing PNIPAM-DNA-DiYO-1
%       particles. Each cell image was classified into one of three
%       categories: single-particle, nulti-particle, or no-particle.
%       Data augmentation (rotation, scaling, reflection) was applied 
%       to improve model robustness.
%
% File outputs:
%       1) Segmented cell images classified into three categories 
%       (Single-particle, Multi-particle, No-particle), stored in 
%       designated folders.
%       2) Trained 3-class ResNet model (.mat file).
%       3) Validation results including accuracy/loss plots and 
%       confusion matrix.
%
% Requirements:
%       MATLAB R2023a or later
%       Image Processing Toolbox
%       ImageJ/Fiji (for initial channel decomposition into TIFF images)
%
% ==========================================================================

clc; clear; warning('off','all'); close all

% Directory (Example directory & Pathes)
main_dir = 'D:/Image Process DNA/LysoView (PAH-AF594)';

Current_folder = fullfile(main_dir);
Sample_folder_info = dir(Current_folder);
Samples = {Sample_folder_info([Sample_folder_info.isdir] & ...
    ~ismember({Sample_folder_info.name}, {'.', '..'})).name};

master_intensity_results = table();
master_null_results      = table();
summary_all_samples      = table(); 

for b = 1:length(Samples)
    % Directory
    SingleCell_path = fullfile(main_dir, Samples{b} ,'Single Cells (BF)');       mkdir(SingleCell_path);
    Green_path      = fullfile(main_dir, Samples{b} ,'Single Cells (Green)');    mkdir(Green_path);
    Red_path        = fullfile(main_dir, Samples{b} ,'Single Cells (Red)');      mkdir(Red_path);
    Cytosol_path    = fullfile(main_dir, Samples{b} ,'Single Cells (Cytosol)');  mkdir(Cytosol_path);

    ImageCount = 0;

    phase_dir  = fullfile(main_dir, Samples{b}, 'Phase');
    blue_dir   = fullfile(main_dir, Samples{b}, 'Blue');
    green_dir  = fullfile(main_dir, Samples{b}, 'Green');
    red_dir    = fullfile(main_dir, Samples{b}, 'Red');

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

%% =================================================================
%  PART 1: Cell Segmentation
% ==================================================================

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
        green_filename = strrep(phase_filename, 'Phase', 'Green');
        red_filename   = strrep(phase_filename, 'Phase', 'Red');
        Img_1 = fullfile(phase_dir, phase_filename);
        Img_2 = fullfile(blue_dir, blue_filename);
        Img_3 = fullfile(green_dir, green_filename);
        Img_4 = fullfile(red_dir, red_filename);

        if ~isfile(Img_1) || ~isfile(Img_2) || ~isfile(Img_3) || ~isfile(Img_4)
            continue;
        end

    %% Reading images

        BF = imread(Img_1);
        Blue = imread(Img_2);
        Green = imread(Img_3);
        Red = imread(Img_4);

        green_full_double = im2double(Green);
        [H_green, W_green] = size(green_full_double);
        border_mask = false(H_green, W_green);
        border_mask(1:10, :) = true;  
        border_mask(end-9:end, :) = true;
        border_mask(:, 1:10) = true; 
        border_mask(:, end-9:end) = true;
        
        border_pixels = green_full_double(border_mask);
        image_background_intensity = median(border_pixels(:));

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

        % Processing Red Images
        I_red_smoothed = imgaussfilt(Red, 3);
        bw_red = imbinarize(I_red_smoothed, 'adaptive',...
            'ForegroundPolarity', 'bright', 'Sensitivity', 0.05);
        bw_clean = imfill(bw_red, 'holes');
        par_mask_clean = imgaussfilt(im2double(bw_clean),5);
        par_mask_clean = par_mask_clean > 0.75;
        D_red = bwdist(~bw_clean);
        D_neg = -D_red;
        markers = imextendedmin(D_neg, 1);
        D_imposed = imimposemin(D_neg, markers);            
        L_red = watershed(D_imposed);             
        par_mask = par_mask_clean;
        par_mask(L_red == 0) = 0;

        label_Mask = bwlabel(mask_clean);
        n_cells = max(label_Mask(:));

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

        label_final = bwlabel(clean_mask);
        image_intensity_results = table();

        for o = 1:max(label_final(:))
            region_final_cell = (label_final == o);

            [row3, col3] = find(label_final == o);
            len3 = max(row3) - min(row3) + 1;
            breadth3 = max(col3) - min(col3) + 1;
            target5 = uint8(zeros([len3 breadth3]));
            target6 = zeros([len3 breadth3]);
            target7 = uint8(zeros([len3 breadth3]));
            target8 = zeros([len3 breadth3]);
            target9 = uint8(zeros([len3 breadth3]));
            sy3 = min(col3) - 1;
            sx3 = min(row3) - 1;

            for p = 1:length(row3)
                x3 = row3(p) - sx3;
                y3 = col3(p) - sy3;
                target5(x3, y3) = BF(row3(p), col3(p)); % Original Bright-field image
                target6(x3, y3) = bw_Nuc_clean(row3(p), col3(p)); % Nuclear Mask image
                target7(x3, y3) = Green(row3(p), col3(p)); % Original Green image
                target8(x3, y3) = par_mask(row3(p), col3(p)); % Red Mask image
                target9(x3, y3) = Red(row3(p), col3(p));
            end

            label_nuc = bwlabel(target6);
            label_red = bwlabel(target8);

            if max(label_nuc(:)) == 1 && max(label_red(:)) == 1
                Stat_redMask = regionprops(label_red, 'Area', 'EquivDiameter', 'Circularity');
                RedMask_area = Stat_redMask(1).Area;
                RedMask_circ = Stat_redMask(1).Circularity;
                tolerance = 0.30;
                Act_area = 450; % ≈ 452 px² (4.8 µm diameter particle ≈72 µm² at 0.2 µm/px)
                min_red_mask_area = round(Act_area*(1-tolerance)); % 315
                max_red_mask_area = round(Act_area*(1+tolerance)); % 585
                min_red_mask_circ = 0.85;

                ImageCount = ImageCount + 1;

                if RedMask_area >= min_red_mask_area && RedMask_area < max_red_mask_area && RedMask_circ >= min_red_mask_circ
                    % imwrite(target5, fullfile(SingleCell_path, sprintf('%d.tif', ImageCount)));
                    % imwrite(target7, fullfile(Green_path, sprintf('%d.tif', ImageCount)));
                    % imwrite(target9, fullfile(Red_path, sprintf('%d.tif', ImageCount)));

                    particle_mask_logical = logical(target8);
                    dilation_radius = round(0.40*(Stat_redMask(1).EquivDiameter / 2));
                    dilation_radius = max(1, dilation_radius);
                    se = strel('disk', dilation_radius);
                    padded_particle_mask = imdilate(particle_mask_logical, se);

                    cell_mask_logical = target7 > 0;

                    cytosol_mask_logical = cell_mask_logical & ~padded_particle_mask;

                    % imwrite(target7 .* uint8(cytosol_mask_logical), fullfile(Cytosol_path, sprintf('%d.tif', ImageCount)));

                    green_img_double = im2double(target7);
                    red_img_double = im2double(target9);

                    Particle_red_intensities   = red_img_double(particle_mask_logical);
                    mean_red_particle   = mean(Particle_red_intensities,'omitnan'); 
                    median_red_particle = median(Particle_red_intensities,'omitnan');
                    std_red_particle    = std(Particle_red_intensities,  'omitnan');

                    Particle_green_intensities = green_img_double(particle_mask_logical);
                    mean_green_particle   = mean(Particle_green_intensities,'omitnan');
                    median_green_particle = median(Particle_green_intensities,'omitnan');
                    std_green_particle    = std(Particle_green_intensities,  'omitnan');

                    Cytosol_green_intensities = green_img_double(cytosol_mask_logical);
                    mean_green_cytosol   = mean(Cytosol_green_intensities,'omitnan');
                    median_green_cytosol = median(Cytosol_green_intensities,'omitnan');
                    std_green_cytosol    = std(Cytosol_green_intensities,  'omitnan');

                    % ==== Pearson correlation in phagosome vs null ====
                    phago_mask = padded_particle_mask & cell_mask_logical;

                    % ===== NEW: Empirical null test (phagosome vs random cytoplasm ROIs) =====
                    N = 1000;  % number of random ROIs
                    true_mean_green = mean(green_img_double(phago_mask),'omitnan');

                    null_means = nan(1,N);
                    phago_area = sum(phago_mask(:));
                    [X,Y] = meshgrid(1:size(green_img_double,2), 1:size(green_img_double,1));
                    [rows, cols] = find(cytosol_mask_logical & ~logical(target6)); % exclude nucleus

                    for i = 1:N
                        idx = randi(length(rows));
                        cx = cols(idx); cy = rows(idx);
                        circle_mask = (X-cx).^2 + (Y-cy).^2 <= sqrt(phago_area/pi)^2;
                        roi_mask = circle_mask & cytosol_mask_logical & ~logical(target6);
                        if sum(roi_mask(:)) >= 5
                            null_means(1,i) = mean(green_img_double(roi_mask),'omitnan');
                        end
                    end

                    percentile_green = 100 * sum(null_means < true_mean_green) / sum(~isnan(null_means));
                    pval_green = (100 - percentile_green)/100;  % one-tailed

                    % Collect per-cell results
                    cell_results = table( ...
                        string(Samples{b}), d, o, RedMask_area, RedMask_circ, ...
                        mean_red_particle, median_red_particle, std_red_particle, ...
                        mean_green_particle, median_green_particle, std_green_particle, ...
                        mean_green_cytosol, median_green_cytosol, std_green_cytosol, ...
                        percentile_green, pval_green, ...
                        'VariableNames', {'Sample','ImageIdx','CellIdx','RedMaskArea','Circularity', ...
                        'MeanRed_Particle','MedianRed_Particle','StdRed_Particle', ...
                        'MeanGreen_Particle','MedianGreen_Particle','StdGreen_Particle', ...
                        'MeanGreen_Cytosol','MedianGreen_Cytosol','StdGreen_Cytosol', ...
                        'Percentile_Green','Pval_Green'} );

                    null_results = table( ...
                        repmat(string(Samples{b}), numel(null_means), 1), ...
                        repmat(d, numel(null_means), 1), ...
                        repmat(o, numel(null_means), 1), ...
                        (1:numel(null_means))', ...
                        null_means', ...
                        'VariableNames', {'Sample','ImageIdx','CellIdx','NullIdx','NullMean'} );

                    cell_results.NullMeans = {null_means};

                    image_intensity_results = [image_intensity_results; cell_results];
                    master_null_results = [master_null_results; null_results];
                end
            end
        end
        master_intensity_results = [master_intensity_results; image_intensity_results];
    end % images loop

    %% =================================================================
    %  SAMPLE-LEVEL SUMMARIES (per sample)
    % ==================================================================
    disp('---------------------------------------------------------');
    fprintf('Calculating sample summaries for: "%s"\n', Samples{b});
    current_sample_name = Samples{b};
    sample_rows_logical = strcmp(master_intensity_results.Sample, current_sample_name);
    current_sample_data = master_intensity_results(sample_rows_logical, :);

    if height(current_sample_data) < 2
        fprintf('  -> Not enough data points (< 2).\n\n');
        continue;
    end

    % --- Print results to Command Window ---
    fprintf('  -> Green enrichment null test: Percentile=%.1f%%, Mean p=%.3g\n', ...
            mean(current_sample_data.Percentile_Green,'omitnan'), ...
            mean(current_sample_data.Pval_Green,'omitnan'));

    summary_row = table( ...Show me the 
        string(current_sample_name), ...
        round(mean(current_sample_data.Percentile_Green,'omitnan'),2), ...
        round(mean(current_sample_data.Pval_Green,'omitnan'),4), ...
        'VariableNames', {'Sample','MeanPercentile','MeanPval'} );

    summary_all_samples = [summary_all_samples; summary_row];
end % samples loop

output_filename = 'Master_Analysis_Results06.xlsx';

if ismember('NullMeans', master_intensity_results.Properties.VariableNames)
    expanded_nulls = cell2mat(master_intensity_results.NullMeans);
    null_colnames = strcat("Null_", string(1:size(expanded_nulls,2)));
    expanded_nulls_table = array2table(expanded_nulls, 'VariableNames', null_colnames);
    master_intensity_results_noNull = removevars(master_intensity_results, 'NullMeans');
    master_intensity_results_expanded = [master_intensity_results_noNull expanded_nulls_table];
    writetable(master_intensity_results_expanded, fullfile(main_dir, output_filename), 'Sheet', 'Measurement_Data');
else
    writetable(master_intensity_results, fullfile(main_dir, output_filename), 'Sheet', 'Measurement_Data');
end

if ~isempty(master_null_results)
    writetable(master_null_results, fullfile(main_dir, output_filename), 'Sheet', 'Null_Distributions');
end

if exist('summary_all_samples','var') && ~isempty(summary_all_samples)
    writetable(summary_all_samples, fullfile(main_dir, output_filename), 'Sheet', 'Summary_Stats');
end
