% 
% clc; close all; clear;
% 
% % --- Read Excel ---
% T = readtable('D:/Image Process DNA/LysoView (PAH-AF594)/Master_Analysis_Results06.xlsx');
% 
% % --- Split samples ---
% s1 = T.Percentile_Green(strcmp(T.Sample,'Sample 1'));
% s2 = T.Percentile_Green(strcmp(T.Sample,'Sample 2'));
% s3 = T.Percentile_Green(strcmp(T.Sample,'Sample 3'));
% 
% % --- Calculate % saturated at 100 ---
% p100 = [mean(s1==100), mean(s2==100), mean(s3==100)] * 100;
% 
% % --- Swarm plot ---
% figure; hold on;
% 
% swarmchart(ones(size(s1)), s1, 100, [0.2 0.6 1], 'filled', ...
%            'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','k');
% swarmchart(2*ones(size(s2)), s2, 100, [1 0.4 0.4], 'filled', ...
%            'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','k');
% swarmchart(3*ones(size(s3)), s3, 100, [0.3 0.8 0.3], 'filled', ...
%            'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','k');
% 
% % --- Formatting ---
% xlim([0.5 3.5]); ylim([96 101]); % zoom in near ceiling
% xticks([1 2 3]); xticklabels({'Sample 1','Sample 2','Sample 3'});
% ylabel('Percentile Green','FontSize',18,'FontWeight','bold');
% % title('Swarm plot of Percentile Green','FontSize',13);
% 
% % --- Annotate % at 100 above each group ---
% for i = 1:3
%     text(i, 100.5, sprintf('%.1f%%\n at 100', p100(i)), ...
%         'HorizontalAlignment','center', ...
%         'FontSize',12,'FontWeight','bold');
% end
% 
% set(gca,'FontSize',18,'LineWidth',1.2,'Box','on');
% grid 'on'
% axis square

%% ============================================================
%  LysoView acidification index (Z-score per phagosome)
%  Robust to "NullMeans" (cell) OR "Null_1..Null_N" (wide) formats
%  Separate swarm and box plots (no notches)
% ============================================================
clc; clear; close all;

% --- Load data ---
xlsx = 'D:/Image Process DNA/LysoView (PAH-AF594)/Master_Analysis_Results06.xlsx';
T = readtable(xlsx);  % if needed: , 'Sheet','Measurement_Data'

% --- Ensure Sample is string for easy comparisons/labels ---
if ~isstring(T.Sample), T.Sample = string(T.Sample); end

% --- Find null layout ---
vnames = T.Properties.VariableNames;
hasNullMeans = ismember('NullMeans', vnames); 
nullWideMask = startsWith(vnames, 'Null_'); 
hasNullWide  = any(nullWideMask);

% --- Sanity: need MeanGreen_Particle to compute z ---
assert(ismember('MeanGreen_Particle', vnames), ...
    'Column "MeanGreen_Particle" not found in your table.');

% --- Compute per-row (per-cell) z-scores ---
Z = nan(height(T),1);

for r = 1:height(T)
    phago_mean = T.MeanGreen_Particle(r);

    if hasNullMeans && ~isempty(T.NullMeans{r})
        % Layout A: a cell containing the 1000 null means
        null_vals = T.NullMeans{r};
    elseif hasNullWide
        % Layout B: expanded columns Null_1..Null_N (same row)
        null_vals = T{r, nullWideMask};
    else
        % No null data present
        null_vals = NaN;
    end

    mu = mean(null_vals, 'omitnan');
    sd = std(null_vals, 0, 'omitnan');
    if sd==0 || isnan(sd)
        Z(r) = NaN;
    else
        Z(r) = (phago_mean - mu) / sd;
    end
end

% --- Assemble z-score table for plotting ---
Tz = table(T.Sample, Z, 'VariableNames', {'Sample','Zscore_Green'});

% --- Colors and labels ---
labs = unique(Tz.Sample, 'stable');
nG = numel(labs);
cols = [0.2 0.6 1; 1 0.4 0.4; 0.3 0.8 0.3; 0.55 0.40 0.80; 0.95 0.7 0.2];
cols = cols(1:nG, :);

%% ============================================================
% 1) Swarm plot (separate figure)
% ============================================================
figure('Color','w','Position',[280 320 700 520]); hold on;

for i = 1:nG
    mask = Tz.Sample == labs(i);
    swarmchart(i*ones(sum(mask),1), Tz.Zscore_Green(mask), ...
        55, cols(i,:), 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.70);
end

yline(0,'--k','Baseline (no enrichment)','LabelHorizontalAlignment','left');
yline(2,':k','z = 2','LabelHorizontalAlignment','left');
yline(3,':k','z = 3','LabelHorizontalAlignment','left');

set(gca,'XTick',1:nG,'XTickLabel',labs,'FontSize',14,'LineWidth',1.2);
ylabel('Z-score (Phagosome LysoView vs cytosol)','FontSize',14,'FontWeight','bold');
title('LysoView acidification index — Swarm','FontSize',14);
grid on; box on; axis square;

% Robust y-limits
allZ = Tz.Zscore_Green(~isnan(Tz.Zscore_Green));
if ~isempty(allZ)
    ylim([min(-0.5, prctile(allZ,2))  max(3.5, prctile(allZ,98))]);
end

%% ============================================================
% 2) Box plot (separate figure, regular — no notches)
% ============================================================
figure('Color','w','Position',[1020 320 700 520]); hold on;

% Boxplot (plain rectangles)
boxplot(Tz.Zscore_Green, Tz.Sample, ...
    'Colors', [0 0 0], ...
    'Symbol', 'o', ...
    'Widths', 0.5, ...
    'Notch', 'off', ...            % <-- regular box, not 'W' / notched
    'BoxStyle', 'outline', ...
    'MedianStyle', 'line');

% Lightly color the boxes per group
hBoxes = findobj(gca,'Tag','Box'); % returned in reverse order
for i = 1:nG
    b = hBoxes(i);
    patch(get(b,'XData'), get(b,'YData'), cols(nG-i+1,:), ...
        'FaceAlpha',0.25, 'EdgeColor','k', 'LineWidth',1.2);
end

yline(0,'--k','Baseline (no enrichment)','LabelHorizontalAlignment','left');
yline(2,':k','z = 2','LabelHorizontalAlignment','left');
yline(3,':k','z = 3','LabelHorizontalAlignment','left');

set(gca,'FontSize',14,'LineWidth',1.2);
ylabel('Z-score (Phagosome LysoView vs cytosol)','FontSize',14,'FontWeight','bold');
title('LysoView acidification index — Box','FontSize',14);
grid on; box on; axis square;

if ~isempty(allZ)
    ylim([min(-0.5, prctile(allZ,2))  max(3.5, prctile(allZ,98))]);
end

%% (Optional) quick console stats
fprintf('\nZ-score summary per sample (mean ± SD | median | %% z>1.96):\n');
for i = 1:nG
    zi = Tz.Zscore_Green(Tz.Sample==labs(i));
    m  = mean(zi,'omitnan'); s = std(zi,'omitnan'); med = median(zi,'omitnan');
    frac2 = 100*mean(zi>1.96,'omitnan');
    fprintf('  %-12s  %.2f ± %.2f  |  med=%.2f  |  %.1f%% > 1.96\n', labs(i), m, s, med, frac2);
end




% clc; close all; clear;
% 
% % Read table
% T = readtable('D:/Image Process DNA/LysoView (PAH-AF594)/Master_Analysis_Results06.xlsx');
% 
% % Split by sample
% s1 = T.Percentile_Green(strcmp(T.Sample,'Sample 1'));
% s2 = T.Percentile_Green(strcmp(T.Sample,'Sample 2'));
% s3 = T.Percentile_Green(strcmp(T.Sample,'Sample 3'));
% 
% %% Histogram overlay (with transparency, normalized)
% figure; hold on;
% histogram(s1,'BinWidth',2,'Normalization','probability',...
%     'FaceColor',[0.2 0.6 1],'FaceAlpha',0.5);
% histogram(s2,'BinWidth',2,'Normalization','probability',...
%     'FaceColor',[1 0.4 0.4],'FaceAlpha',0.5);
% histogram(s3,'BinWidth',2,'Normalization','probability',...
%     'FaceColor',[0.3 0.8 0.3],'FaceAlpha',0.5);
% xlabel('Percentile Green'); ylabel('Probability');
% legend({'Sample 1','Sample 2','Sample 3'});
% title('Normalized Histograms');
% xlim([90 101]); % zoom into ceiling region
% grid on;
% 
% %% ECDF overlay with styled lines
% figure; hold on;
% 
% [f1,x1] = ecdf(s1);
% [f2,x2] = ecdf(s2);
% [f3,x3] = ecdf(s3);
% 
% plot(x1,f1,'-','LineWidth',2,'Color',[0.2 0.6 1]);
% plot(x2,f2,'-','LineWidth',2,'Color',[1 0.4 0.4]);
% plot(x3,f3,'-','LineWidth',2,'Color',[0.3 0.8 0.3]);
% 
% legend({'Sample 1','Sample 2','Sample 3'},'Location','southeast');
% xlabel('Percentile Green'); ylabel('Cumulative Probability');
% title('ECDF Comparison');
% xlim([90 101]); % zoom into ceiling region
% grid on;
% 
% % Example categorization
% edges = [0 95 99.9 100.1]; % bins
% cats1 = discretize(s1,edges);
% cats2 = discretize(s2,edges);
% cats3 = discretize(s3,edges);
% 
% counts = [histcounts(s1,edges,'Normalization','probability'); ...
%           histcounts(s2,edges,'Normalization','probability'); ...
%           histcounts(s3,edges,'Normalization','probability')];
% 
% figure;
% bar(counts,'stacked');
% xticklabels({'Sample 1','Sample 2','Sample 3'});
% legend({'<95','95–99.9','=100'},'Location','northoutside','Orientation','horizontal');
% ylabel('Proportion');
% title('Distribution categories');
% 
% b = bar(counts,'stacked');
% xticklabels({'Sample 1','Sample 2','Sample 3'});
% 
% figure;
% 
% % --- Left: Stacked proportion bar plot ---
% subplot(1,2,1);
% bar(counts,'stacked');
% xticklabels({'Sample 1','Sample 2','Sample 3'});
% ylabel('Proportion');
% legend({'<95','95–99.9','=100'},'Location','northoutside','Orientation','horizontal');
% title('Distribution categories');
% 
% % Add percentage labels
% for i = 1:3
%     for j = 1:3
%         if counts(i,j) > 0
%             text(i, sum(counts(i,1:j))-counts(i,j)/2, ...
%                 sprintf('%.1f%%', counts(i,j)*100), ...
%                 'HorizontalAlignment','center','FontSize',8,'Color','k');
%         end
%     end
% end
% 
% % --- Right: Jittered dot plot ---
% subplot(1,2,2); hold on;
% scatter(ones(size(s1)),s1,15,[0.2 0.6 1],'filled','jitter','on','jitterAmount',0.15);
% scatter(2*ones(size(s2)),s2,15,[1 0.4 0.4],'filled','jitter','on','jitterAmount',0.15);
% scatter(3*ones(size(s3)),s3,15,[0.3 0.8 0.3],'filled','jitter','on','jitterAmount',0.15);
% xlim([0.5 3.5]); ylim([90 101]); % zoom in to show differences near ceiling
% xticks([1 2 3]); xticklabels({'Sample 1','Sample 2','Sample 3'});
% ylabel('Percentile Green');
% title('Jittered dot plot (zoomed)');

% clc; close all; clear;
% 
% % Read Excel
% T = readtable('D:/Image Process DNA/LysoView (PAH-AF594)/Master_Analysis_Results06.xlsx');
% 
% % Extract groups
% s1 = T.Percentile_Green(strcmp(T.Sample,'Sample 1'));
% s2 = T.Percentile_Green(strcmp(T.Sample,'Sample 2'));
% s3 = T.Percentile_Green(strcmp(T.Sample,'Sample 3'));
% 
% % Calculate % at 100
% p100 = [mean(s1==100), mean(s2==100), mean(s3==100)] * 100;
% 
% figure;
% 
% % --- Left: Lollipop plot ---
% subplot(1,2,1); hold on;
% x = 1:3;
% for i = 1:3
%     plot([x(i) x(i)],[0 p100(i)],'k-','LineWidth',1.5); % stick
%     plot(x(i),p100(i),'o','MarkerSize',8,'MarkerFaceColor',[0 0.45 0.74],'MarkerEdgeColor','k'); % dot
% end
% xlim([0.5 3.5]); ylim([0 105]);
% xticks(1:3); xticklabels({'Sample 1','Sample 2','Sample 3'});
% ylabel('% at 100%');
% title('Proportion saturated at 100%');
% 
% % --- Right: Jittered dot plot (zoomed) ---
% subplot(1,2,2); hold on;
% scatter(ones(size(s1)),s1,15,[0.2 0.6 1],'filled','jitter','on','jitterAmount',0.15);
% scatter(2*ones(size(s2)),s2,15,[1 0.4 0.4],'filled','jitter','on','jitterAmount',0.15);
% scatter(3*ones(size(s3)),s3,15,[0.3 0.8 0.3],'filled','jitter','on','jitterAmount',0.15);
% xlim([0.5 3.5]); ylim([90 101]);
% xticks([1 2 3]); xticklabels({'Sample 1','Sample 2','Sample 3'});
% ylabel('Percentile Green');
% title('Non-saturated values (zoomed)');
% 
% figure; hold on;
% scatter(ones(size(s1)),s1,15,[0.2 0.6 1],'filled','jitter','on','jitterAmount',0.15);
% scatter(2*ones(size(s2)),s2,15,[1 0.4 0.4],'filled','jitter','on','jitterAmount',0.15);
% scatter(3*ones(size(s3)),s3,15,[0.3 0.8 0.3],'filled','jitter','on','jitterAmount',0.15);
% 
% xlim([0.5 3.5]); ylim([90 101]);
% xticks([1 2 3]); xticklabels({'Sample 1','Sample 2','Sample 3'});
% ylabel('Percentile Green');
% title('Raw values with % saturated');
% 
% % Annotate % at 100 above groups
% for i = 1:3
%     text(i,101,sprintf('%.1f%% at 100',p100(i)),...
%         'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
% end
% 
% figure; hold on;
% swarmchart(ones(size(s1)),s1,15,[0.2 0.6 1],'filled');
% swarmchart(2*ones(size(s2)),s2,15,[1 0.4 0.4],'filled');
% swarmchart(3*ones(size(s3)),s3,15,[0.3 0.8 0.3],'filled');
% xlim([0.5 3.5]); ylim([90 101]);
% xticks([1 2 3]); xticklabels({'Sample 1','Sample 2','Sample 3'});
% ylabel('Percentile Green');
% title('Swarm plot of values');