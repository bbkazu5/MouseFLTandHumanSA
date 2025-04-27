% Extract flight experiment data
load('C:\Users\khan332\Documents\GitHub\MouseFLTandHumanSA\BehaviourData\GaitData.mat')
flightData = GaitData(GaitData.Treatment == 'Flight', :);
groundData = GaitData(GaitData.Treatment == 'GC', :);
%%
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

% % Grab ground control data pre- and post-
% GC_pre = flightData(groundData.Time == 'Pre', :);
% GC_pre = flight_pre(GC_pre.Limb1 == 'LH', :);
% GC_post = flightData(groundData.Time == 'Post', :);
% GC_post = flight_post(GC_post.Limb1 == 'LH', :);


effect_sizes = table();
for v = 1:length(vars)
    [~,~,~,stats] = ttest2(flight_pre{:,vars{v}}, flight_post{:,vars{v}});
    d = abs(mean(flight_post{:,vars{v}}) - mean(flight_pre{:,vars{v}})) / std([flight_pre{:,vars{v}}; flight_post{:,vars{v}}]);
    effect_sizes.(vars{v}) = d;
end

% Convert to sorted vector
[~,idx] = sort(effect_sizes.Variables, 'descend');
sorted_vars = vars(idx);
effectData = cell2mat(table2cell(effect_sizes));
effectDataSorted = effectData(idx);

%% Plot Cohens D
figure('Position', [100, 100, 800, 800])
barh(effectDataSorted)  % Show top 20 only
xlabel("Cohen's d Effect Size")
title('Behavioral Parameter Changes: Flight vs Pre-Post')
% Remove common prefixes/suffixes
% Suppose sorted_vars is your cell array of variable names
cleanLabels = sorted_vars;

% Replace underscores with spaces
cleanLabels = strrep(cleanLabels, '_', ' ');

% (Optional) Remove units like 'centimeters', 'seconds', 'degrees', 'percent', 'ratio', 'coefficient of variance', etc.
units = {'centimeters','seconds','degrees','percent','ratio','coefficient of variance','squared'};
for i = 1:numel(units)
    cleanLabels = regexprep(cleanLabels, [' ' units{i}], '');
end

% (Optional) Capitalize first letter of each word
for i = 1:numel(cleanLabels)
    words = split(cleanLabels{i});
    words = cellfun(@(w) [upper(w(1)), lower(w(2:end))], words, 'UniformOutput', false);
    cleanLabels{i} = strjoin(words, ' ');
end
set(gca, 'YTick', 1:length(cleanLabels), 'YTickLabel', cleanLabels)
%% Plot Clustered data
% Calculate correlations
corr_matrix = corr([flight_pre{:,sorted_vars}; flight_post{:,sorted_vars}]);

cg = clustergram(corr_matrix, ...
    'Colormap', parula, ...
    'RowLabels', cleanLabels, ...
    'ColumnLabels', cleanLabels, ...
    'Standardize', 'none', ...
    'RowPdist', 'correlation', ...
    'ColumnPdist', 'correlation');
cg.XDisplayLabels = cleanLabels;
cg.YDisplayLabels = cleanLabels;
cg.Title = 'Behavioral Parameter Correlations';
%% Plot most prominant behaviours
Flight = struct();
GC = struct();
[Flight.preStrideLength,Flight.postStrideLength] = getStride_Length(flightData);
[GC.preStrideLength,GC.postStrideLength] = getStride_Length(groundData);


% Stride Variability
[Flight.preStrideVar,Flight.postStrideVar] = getStride_SwingTime(flightData);
[GC.preStrideVar,GC.postStrideVar] = getStride_SwingTime(groundData);


% Perform paired t-test
[h, p, ci, stats] = ttest(Flight.preStrideLength, Flight.postStrideLength);
fprintf('Flight Paired t-test:\n p = %.4f, t(%d) = %.2f\n', p, stats.df, stats.tstat);
[h, p, ci, stats] = ttest(GC.preStrideLength, GC.postStrideLength);
fprintf('Flight Paired t-test:\n p = %.4f, t(%d) = %.2f\n', p, stats.df, stats.tstat);
%%

labels = {'Pre-Flight', 'Post-Flight', 'Pre-GC','Post-GC'};
data = [Flight.preStrideLength,Flight.postStrideLength,GC.preStrideLength(1:17),GC.postStrideLength(1:17)];
figure,
makeBoxplot(data, labels)

labels = {'Pre-Flight', 'Post-Flight', 'Pre-GC','Post-GC'};
data = [Flight.preStrideVar,Flight.postStrideVar,GC.preStrideVar(1:17),GC.postStrideVar(1:17)];
figure,
makeBoxplot(data, labels)
ylabel('Swing Time (s)')
%% LOCAL FUNCTIONS
function [preStride,postStride] = getStride_Length(data)
% Separate pre- and post-flight measurements
preData = data(data.Time == 'Pre', :);
preData = preData(preData.Limb == 'Left Hind', :);
postData = data(data.Time == 'Post', :);
postData = postData(postData.Limb == 'Left Hind', :);

% Get paired stride lengths (assuming same subjects in both groups)
preStride = preData.StrideLength_centimeters(:);
postStride = postData.StrideLength_centimeters(:);
assert(size(preStride,1)==size(postStride,1))
end

function [preStride,postStride] = getStride_SwingTime(data)
% Separate pre- and post-flight measurements
preData = data(data.Time == 'Pre', :);
preData = preData(preData.Limb == 'Left Hind', :);
postData = data(data.Time == 'Post', :);
postData = postData(postData.Limb == 'Left Hind', :);

% Get paired stride lengths (assuming same subjects in both groups)
preStride = preData.SwingTime_seconds(:);
postStride = postData.SwingTime_seconds(:);
assert(size(preStride,1)==size(postStride,1))
end

function makeBoxplot(data, groupLabels)
% makeBoxplot(data, groupLabels)
% data: matrix, each column is a group
% groupLabels: cell array of strings, one per column

if nargin < 2
    groupLabels = strcat("Group ", string(1:size(data,2)));
end

% Generate distinguishable colors
boxColors = lines(size(data,2));

% Boxplot with custom colors
h = boxplot(data, ...
    'Labels', groupLabels, ...
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
set(gca, 'FontSize', 13, 'LineWidth', 1.2);
title('Spaceflight Effects on Mouse Gait', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gcf, 'Color', 'w');

% --- ANOVA ---
[p, tbl, stats] = anova1(data, groupLabels, 'off'); % 'off' suppresses the default ANOVA plot

% Extract F-statistic from ANOVA table
F_value = tbl{2,5}; % Row 2 (Treatment), Column 5 (F-statistic)



% Update fprintf statement
fprintf('One-way ANOVA:\nF(%d,%d) = %.2f, p = %.4f\n',...
        tbl{2,3}, tbl{3,3}, F_value, p);

% Annotate plot with p-value
yLimits = ylim;
xPosition = mean(xlim); % Centered
yPosition = yLimits(2) - 0.07*diff(yLimits); % 7% below top
% Create formatted annotation string
annotationText = sprintf('One-way ANOVA:\nF(%d,%d) = %.2f\np = %.4f',...
                        tbl{2,3}, tbl{3,3}, F_value, p);
text(xPosition, yPosition, annotationText, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'BackgroundColor', [1 1 1 0.7], ...
    'Margin', 3);
end
