flightData = GaitData(GaitData.Treatment == 'Flight', :);
groundData = GaitData(GaitData.Treatment == 'GC', :);

% List of columns to exclude (non-behavioral/metadata)
exclude = {'SourceName','SampleName','RFID','Limb','Limb1','Time','Treatment'};

% Get all variable names
allVars = GaitData.Properties.VariableNames;

% Remove excluded variables
vars = setdiff(allVars, exclude);

% Grab preflight and post flight data
flight_pre = flightData(flightData.Time == 'Pre', :);
flight_pre = flight_pre(flight_pre.Limb1 == 'LH', :);
flight_post = flightData(flightData.Time == 'Post', :);
flight_post = flight_post(flight_post.Limb1 == 'LH', :);

% Grab ground control data pre- and post- (fixed these lines)
GC_pre = groundData(groundData.Time == 'Pre', :);
GC_pre = GC_pre(GC_pre.Limb1 == 'LH', :);
GC_post = groundData(groundData.Time == 'Post', :);
GC_post = GC_post(GC_post.Limb1 == 'LH', :);

%% Combine data for PCA
% Create group labels
flight_pre_labels = repmat({'Flight_Pre'}, size(flight_pre, 1), 1);
flight_post_labels = repmat({'Flight_Post'}, size(flight_post, 1), 1);
GC_pre_labels = repmat({'GC_Pre'}, size(GC_pre, 1), 1);
GC_post_labels = repmat({'GC_Post'}, size(GC_post, 1), 1);

% Combine data
all_data = [flight_pre; flight_post; GC_pre; GC_post];
group_labels = [flight_pre_labels; flight_post_labels; GC_pre_labels; GC_post_labels];

% Extract behavioral variables
X = all_data{:, vars};

% Standardize data (mean=0, std=1)
X_std = zscore(X);

%% Perform PCA
[coeff, score, latent, ~, explained] = pca(X_std);

%% Create PCA visualization
figure('Position', [100, 100, 1000, 800], 'Color', 'w');

% Create a 2x2 subplot layout
subplot(2, 2, 1);
% Plot PC1 vs PC2 scatter
gscatter(score(:,1), score(:,2), group_labels, 'rgbm', '.', 20);
xlabel('PC1 (Variability explained: ' + string(round(explained(1),1)) + '%)');
ylabel('PC2 (Variability explained: ' + string(round(explained(2),1)) + '%)');
title('PCA of Gait Parameters (PC1 vs PC2)');
grid on;
axis square
% Add 95% confidence ellipses for each group
hold on;
groups = unique(group_labels);
colors = 'rgbm';
for i = 1:length(groups)
    idx = strcmp(group_labels, groups{i});
    h = error_ellipse(cov(score(idx,1:2)), mean(score(idx,1:2)), colors(i));
    set(h, 'LineWidth', 2);
end
hold off;
legend('Location', 'eastoutside');
axis square
subplot(2, 2, 2);
% Plot PC1 vs PC3 scatter
gscatter(score(:,1), score(:,3), group_labels,'rgbm', '.', 20);
xlabel('PC1 (Variability explained: ' + string(round(explained(1),1)) + '%)');
ylabel('PC3 (Variability explained: ' + string(round(explained(3),1)) + '%)');
title('PCA of Gait Parameters (PC1 vs PC3)');
grid on;
legend('Location', 'eastoutside');
axis square
% Scree plot
subplot(2, 2, 3);

bar(explained(1:min(10,length(explained))));hold on

plot(cumsum(explained(1:min(10,length(explained)))),'ko-','Linewidth',2)
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('PCA Variance Explained');
grid on;
axis square
% Variable contributions (loadings) to PC1 and PC2
subplot(2, 2, 4);
% Select top contributing variables for readability
[~, idx] = sort(abs(coeff(:,1)) + abs(coeff(:,2)), 'descend');
top_vars = min(4, length(vars)); % Show top 8 variables or fewer
biplot(coeff(idx(1:top_vars),1:2), 'Scores', score(:,1:2), 'VarLabels', vars(idx(1:top_vars)));
title('Variable Contributions to PC1 and PC2');
grid on;
axis square
% Create a colormap for the groups
colormap = [1 0 0; 0 1 0; 0 0 1; 1 0 1];
h = findobj(gca, 'Type', 'Line');
for i = 1:length(h)
    if strcmp(get(h(i), 'Marker'), 'o') || strcmp(get(h(i), 'Marker'), '+') || ...
       strcmp(get(h(i), 'Marker'), '*') || strcmp(get(h(i), 'Marker'), 'x')
        group_idx = strcmp(get(h(i), 'DisplayName'), groups);
        if any(group_idx)
            set(h(i), 'Color', colormap(find(group_idx),:));
        end
    end
end

%% Statistical test for group differences in PC space
% MANOVA on first 3 PCs
X_manova = score(:,1:3);
[d,p,stats] = manova1(X_manova, group_labels);
disp(['MANOVA p-value: ' num2str(p')]);

% Pairwise Hotelling's T-squared tests
fprintf('\nPairwise Hotelling T-squared tests:\n');
group_pairs = nchoosek(1:length(groups), 2);
for i = 1:size(group_pairs, 1)
    g1 = group_pairs(i,1);
    g2 = group_pairs(i,2);
    idx1 = strcmp(group_labels, groups{g1});
    idx2 = strcmp(group_labels, groups{g2});
    [~, pval, ~] = HotellingT2(score(idx1,1:3), score(idx2,1:3));
    fprintf('%s vs %s: p = %.4f\n', groups{g1}, groups{g2}, pval);
end
%% Follow-up: Compare specific behavioral parameters identified by PCA
% Identify top contributing variables to the PCs that separate Flight_Post
[~, topPC1vars] = sort(abs(coeff(:,1)), 'descend');
[~, topPC2vars] = sort(abs(coeff(:,2)), 'descend');

% Select top 3 variables from each PC
top_vars_idx = unique([topPC1vars(1:3); topPC2vars(1:3)]);
top_vars_names = vars(top_vars_idx);

% Create boxplots for these variables
figure('Position', [100, 100, 1200, 600]);
for i = 1:length(top_vars_idx)
    subplot(2, ceil(length(top_vars_idx)/2), i);
    
    % Extract the variable data for each group
    var_idx = top_vars_idx(i);
    var_name = vars{var_idx};
    
    % Combine data and create group labels
    var_data = [flight_pre{:,var_name}; flight_post{:,var_name}; ...
                GC_pre{:,var_name}; GC_post{:,var_name}];
    
    % Create boxplot
    makeBoxplot(var_data, group_labels)
    ylabel(strrep(var_name, '_', ' '));
    title(['Distribution of ' strrep(var_name, '_', ' ')]);
    
    % Run ANOVA for this variable
    [p, tbl, stats] = anova1(var_data, group_labels, 'off');
    
    % Display p-value on the plot
    if p < 0.001
        p_text = 'p < 0.001';
    else
        p_text = ['p = ' num2str(round(p, 3))];
    end
    text(0.1, 0.9, p_text, 'Units', 'normalized');
    
    % Run post-hoc test if ANOVA is significant
    if p < 0.05
        c = multcompare(stats, 'Display', 'off');
        % Check if Flight_Post is different from Flight_Pre
        idx_flight = find(strcmp('Flight_Post', groups)) - find(strcmp('Flight_Pre', groups));
        if any(c(:,1) == 1 & c(:,2) == 1+idx_flight & c(:,6) < 0.05)
            text(0.1, 0.8, 'Flight Post ≠ Flight Pre', 'Units', 'normalized', 'Color', 'r');
        end
    end
end

%% Combined PCA and Behavioral Variable Boxplots Figure

% Prepare figure with tiledlayout
figure('Position', [100, 100, 700, 900], 'Color', 'w');
tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: PC1 vs PC2 scatter
nexttile([1 2]);
gscatter(score(:,1), score(:,2), group_labels, 'rgbm', '.', 20);
xlabel(['PC1 (Variability explained: ' num2str(round(explained(1),1)) '%)']);
ylabel(['PC2 (Variability explained: ' num2str(round(explained(2),1)) '%)']);
title('PCA of Gait Parameters (PC1 vs PC2)');
grid on;
axis square;
hold on;
groups = unique(group_labels);
colors = 'grmb';
for i = 1:length(groups)
    idx = strcmp(group_labels, groups{i});
    h = error_ellipse(cov(score(idx,1:2)), mean(score(idx,1:2)), colors(i));
    set(h, 'LineWidth', 2);
end
hold off;
%legend('Location', 'eastoutside');
legend off

% Panel 2: PC1 vs PC3 scatter
nexttile([1 2]);
gscatter(score(:,1), score(:,3), group_labels, 'rgbm', '.', 20);
xlabel(['PC1 (Variability explained: ' num2str(round(explained(1),1)) '%)']);
ylabel(['PC3 (Variability explained: ' num2str(round(explained(3),1)) '%)']);
title('PCA of Gait Parameters (PC1 vs PC3)');
grid on;
hold on;
groups = unique(group_labels);
colors = 'grmb';
for i = 1:length(groups)
    idx = strcmp(group_labels, groups{i});
    h = error_ellipse(cov(score(idx,[1,3])), mean(score(idx,[1,3])), colors(i));
    set(h, 'LineWidth', 2);
end
hold off;
%legend('Location', 'eastoutside');
axis square;
legend off

% Panel 3: Scree plot
nexttile([1 2]);
bar(explained(1:min(10,length(explained)))); hold on;
plot(cumsum(explained(1:min(10,length(explained)))), 'ko-', 'LineWidth', 2);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('PCA Variance Explained');
grid on;
axis square;

% Panel 4: Variable contributions biplot
nexttile([1 2]);
[~, idx] = sort(abs(coeff(:,1)) + abs(coeff(:,2)), 'descend');
top_vars = min(4, length(vars));
biplot(coeff(idx(1:top_vars),1:2), 'Scores', score(:,1:2), 'VarLabels', vars(idx(1:top_vars)));
title('Variable Contributions to PC1 and PC2');
grid on;
axis square;
box on


% Identify top contributing variables to the PCs
[~, topPC1vars] = sort(abs(coeff(:,1)), 'descend');
[~, topPC2vars] = sort(abs(coeff(:,2)), 'descend');

% Select top 3 variables from each PC (up to 6 unique variables)
top_vars_idx = unique([topPC1vars(1:5); topPC2vars(1:5);]);
top_vars_idx = top_vars_idx([4,7,8,10]);
top_vars_names = vars(top_vars_idx);

% Create boxplots for these variables below
num_boxplots = length(top_vars_idx);
rows = ceil(num_boxplots / 4);
cols = min(4, num_boxplots);

for i = 1:num_boxplots
    % Calculate position in the grid
    row = 1 + ceil(i / cols);
    col = mod(i-1, cols) + 1;
    
    % Create tile for this boxplot
    nexttile(6*row + col - 4);
    
    % Extract the variable data for each group
    var_idx = top_vars_idx(i);
    var_name = vars{var_idx};
    
    % Combine data and create group labels
    var_data = [flight_pre{:,var_name}; flight_post{:,var_name}; ...
                GC_pre{:,var_name}; GC_post{:,var_name}];
    
    % Create boxplot
    makeBoxplot(var_data, group_labels);
    ylabel(strrep(var_name, '_', ' '));
    title(['Distribution of ' strrep(var_name, '_', ' ')]);
    
    % Run ANOVA for this variable
    [p, tbl, stats] = anova1(var_data, group_labels, 'off');
    
    % Display p-value on the plot
    if p < 0.001
        p_text = 'p < 0.001';
    else
        p_text = ['p = ' num2str(round(p, 3))];
    end
    text(0.1, 0.9, p_text, 'Units', 'normalized');
    
    % Run post-hoc test if ANOVA is significant
    if p < 0.05
        c = multcompare(stats, 'Display', 'off');
        % Check if Flight_Post is different from Flight_Pre
        idx_flight = find(strcmp('Flight_Post', groups)) - find(strcmp('Flight_Pre', groups));
        if any(c(:,1) == 1 & c(:,2) == 1+idx_flight & c(:,6) < 0.05)
            text(0.1, 0.8, 'Flight Post ≠ Flight Pre', 'Units', 'normalized', 'Color', 'r');
        end
    end

end

% Add overall title
title(tl, 'PCA of Gait Parameters and Top Contributing Variables', 'FontSize', 14);

% Add summary statistics as annotation
annotation('textbox', [0.01, 0.01, 0.3, 0.03], 'String', ...
    ['MANOVA p-value: ' num2str(p)], ...
    'EdgeColor', 'none', 'FontSize', 9);

%% Helper function for confidence ellipses
function h = error_ellipse(C, mu, color)
    [eigvec, eigval] = eig(C);
    [~, ord] = sort(diag(eigval), 'descend');
    eigvec = eigvec(:, ord);
    eigval = diag(eigval);
    eigval = eigval(ord);
    
    % Calculate the 95% confidence interval
    chisquare_val = 5.991; % 95% confidence for 2 degrees of freedom
    theta = linspace(0, 2*pi, 100);
    a = sqrt(chisquare_val * eigval(1));
    b = sqrt(chisquare_val * eigval(2));
    
    % Parametric equation of the ellipse
    ellipse = [a*cos(theta); b*sin(theta)];
    % Rotate and translate the ellipse to the actual center
    ellipse = eigvec * ellipse + mu';
    
    % Plot the ellipse
    h = plot(ellipse(1,:), ellipse(2,:), color, 'LineWidth', 2);
end

%% Helper function for Hotelling's T-squared test
function [T2, p, stats] = HotellingT2(X1, X2)
    n1 = size(X1, 1);
    n2 = size(X2, 1);
    p = size(X1, 2);
    
    % Calculate means
    mu1 = mean(X1);
    mu2 = mean(X2);
    
    % Calculate pooled covariance
    S1 = cov(X1);
    S2 = cov(X2);
    Sp = ((n1-1)*S1 + (n2-1)*S2) / (n1 + n2 - 2);
    
    % Calculate T^2 statistic
    T2 = (n1*n2)/(n1+n2) * (mu1-mu2) / Sp * (mu1-mu2)';
    
    % Calculate F statistic
    F = (n1+n2-p-1)/((n1+n2-2)*p) * T2;
    
    % p-value
    p = 1 - fcdf(F, p, n1+n2-p-1);
    
    stats.T2 = T2;
    stats.F = F;
    stats.df1 = p;
    stats.df2 = n1+n2-p-1;
end

function makeBoxplot(data, groupLabels)
% makeBoxplot(data, groupLabels)
% data: matrix, each column is a group
% groupLabels: cell array of strings, one per column

if nargin < 2
    groupLabels = strcat("Group ", string(1:size(data,2)));
end

% Generate distinguishable colors
boxColors = lines(length(unique(groupLabels)));

% Boxplot with custom colors
h = boxplot(data, groupLabels, ...
    'Colors', boxColors, ...
    'BoxStyle', 'filled', ...
    'MedianStyle', 'line', ...
    'Symbol', 'k+', ...
    'Widths', 0.6, ...
    'PlotStyle', 'compact');

% Set box colors
boxes = findobj(gca, 'Tag', 'Box');
for j = 1:length(boxes)
    patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), boxColors(j,:), ...
        'FaceAlpha', 0.6, 'EdgeColor', 'none');
end

% Beautify plot
set(findobj(gca,'Tag','Median'),'Color','k','LineWidth',2);
set(findobj(gca,'Tag','Whisker'),'LineWidth',1.5);
set(findobj(gca,'Tag','Upper Whisker'),'LineWidth',1.5);
set(findobj(gca,'Tag','Lower Whisker'),'LineWidth',1.5);
set(findobj(gca,'Tag','Upper Adjacent Value'),'LineWidth',1.5);
set(findobj(gca,'Tag','Lower Adjacent Value'),'LineWidth',1.5);
ylabel('Stride Length (cm)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 7, 'LineWidth', 1.2);
title('Spaceflight Effects on Mouse Gait', 'FontSize', 8, 'FontWeight', 'bold');
grid on;
set(gcf, 'Color', 'w');
end