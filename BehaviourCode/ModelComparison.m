
%%
clear
CPGModel
close all
plotCPGModel
healthy_data.flexor_neurons = flexor_smooth;
healthy_data.extensor_neurons = extensor_smooth;
healthy_data.swing_times = swing_times;
healthy_data.swing_stance_ratio = normrnd(swing_stance_ratio,swing_stance_ratio/10, 15,1)-20;
healthy_data.stance_times = stance_times;
healthy_data.phase_diff = phase_diff;
healthy_data.stride_lengths = stride_lengths;

CPGModel_Sarcopenia
close all

plotCPGModel
sarcopenia_data.flexor_neurons = flexor_smooth;
sarcopenia_data.extensor_neurons = extensor_smooth;
sarcopenia_data.swing_times = swing_times;
sarcopenia_data.swing_stance_ratio = normrnd(swing_stance_ratio,swing_stance_ratio/10, 15,1);
sarcopenia_data.stance_times = stance_times;
sarcopenia_data.phase_diff = phase_diff;
sarcopenia_data.stride_lengths = stride_lengths;
%% Plot and compare both

createGaitComparisonFigure(healthy_data, sarcopenia_data)

%%% LOCAL FUNCTIONS
function createGaitComparisonFigure(healthy_data, sarcopenia_data)
    % Create figure with 2x3 subplot layout
    figure('Position', [100 100 1200 800], 'Color', 'w');
    
    % Color scheme
    healthy_color = [0.4 0.5 0.8];  % Blue tones
    sarcopenia_color = [0.8 0.5 0.4];  % Red-orange tones
    
    %% Panel 1: Swing Time Distribution Comparison
    subplot(2,3,1);
    histogram(healthy_data.swing_times, 'BinWidth', 0.2, 'Normalization', 'probability',...
        'FaceColor', healthy_color, 'EdgeColor', 'none');
    hold on;
    histogram(sarcopenia_data.swing_times, 'BinWidth', 0.2, 'Normalization', 'probability',...
        'FaceColor', sarcopenia_color, 'EdgeColor', 'none');
    xlabel('Swing Time (ms)');
    ylabel('Probability');
    title('A: Swing Time Distribution');
    legend({'Healthy', 'Sarcopenia'}, 'Location', 'northwest');
    box off;
    
    %% Panel 2: Stride Length Distribution Comparison
    subplot(2,3,2);
    histogram(healthy_data.stride_lengths, 'BinWidth', 50, 'Normalization', 'probability',...
        'FaceColor', healthy_color, 'EdgeColor', 'none');
    hold on;
    histogram(sarcopenia_data.stride_lengths, 'BinWidth', 50, 'Normalization', 'probability',...
        'FaceColor', sarcopenia_color, 'EdgeColor', 'none');
    xlabel('Stride Length (ms)');
    ylabel('Probability Density');
    title('B: Stride Length Distribution');
    box off;
    
    %% Panel 3: Swing-Stance Ratio Comparison
    subplot(2,3,3);
    groupLabels = {'Healthy', 'Sarcopenia'};
    makeBoxplot([healthy_data.swing_stance_ratio, sarcopenia_data.swing_stance_ratio], groupLabels)
    hold off;
    ylabel('Swing/Stance Ratio');
    title('C: Swing-Stance Ratio Over Cycles');
    box off;
    
    %% Panel 4: Swing vs Stride Relationship
    subplot(2,3,4);
    hold on;
    plotSwingStrideRelationship(healthy_data, healthy_color);
    plotSwingStrideRelationship(sarcopenia_data, sarcopenia_color);
    xlabel('Swing Time (ms)');
    ylabel('Stride Length (ms)');
    title('D: Swing Time vs Stride Length');
    legend({'Healthy', 'Sarcopenia'}, 'Location', 'northwest');
    box off;
    
    %% Panel 5: Phase Difference Distribution
    subplot(2,3,5);
    polarhistogram(healthy_data.phase_diff, 20, 'FaceColor', healthy_color,...
        'EdgeColor', 'none', 'Normalization', 'probability');
    hold on;
    polarhistogram(sarcopenia_data.phase_diff, 20, 'FaceColor', sarcopenia_color,...
        'EdgeColor', 'none', 'Normalization', 'probability');
    title('E: Phase Difference Distribution');
    
    %% Panel 6: Neural Activity Comparison (Raster Plot Example)
    subplot(2,3,6);
    % Example neural activity plot - modify based on your actual CPG data
    plotNeuralActivity(healthy_data.flexor_neurons, healthy_data.extensor_neurons, healthy_color);
    hold on;
    plotNeuralActivity(sarcopenia_data.flexor_neurons, sarcopenia_data.extensor_neurons, sarcopenia_color);
    xlabel('Time (ms)');
    ylabel('Neuron ID');
    title('F: CPG Neural Activity');
    box off;
    
    % Add overall annotation
    annotation('textbox', [0.1 0.95 0.8 0.05], 'String',...
        'Gait Parameter Comparison: Healthy vs Sarcopenia Conditions',...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
end

function plotSwingStrideRelationship(data, color)
    min_length = min(length(data.stance_times), length(data.swing_times));
    stride_from_phases = data.swing_times(1:min_length) + data.stance_times(1:min_length);
    
    scatter(data.swing_times(1:min_length), stride_from_phases, 40, color,...
        'filled', 'MarkerFaceAlpha', 0.7);
    
    % Add regression line
    p = polyfit(data.swing_times(1:min_length), stride_from_phases, 1);
    x_range = linspace(min(data.swing_times), max(data.swing_times), 100);
    y_fit = polyval(p, x_range);
    plot(x_range, y_fit, 'Color', color, 'LineWidth', 2);
    
    % Add correlation coefficient
    R = corrcoef(data.swing_times(1:min_length), stride_from_phases);
    text(0.1, 0.9, ['r = ', num2str(round(R(1,2),2))],...
        'Units', 'normalized', 'Color', color);
end

function plotNeuralActivity(flexors, extensors, color)
    % Example neural activity plotting - modify based on your data structure
    time = 1:size(flexors,2);
    plot(time, extensors, '--','Color', color, 'LineWidth', 2);
    plot(time, flexors, 'Color', color, 'LineWidth', 2);
    xlim([0 100])
end

function makeBoxplot(data, groupLabels)
% makeBoxplot(data, groupLabels)
% data: matrix, each column is a group
% groupLabels: cell array of strings, one per column

if nargin < 2
    groupLabels = strcat("Group ", string(1:size(data,2)));
end

% Generate distinguishable colors
boxColors = lines(size(data, 2));

% Create basic boxplot
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

% Add raw data points with jitter
hold on;
for group = 1:size(data, 2)
    % Get x-position for this group
    xPos = group;
    
    % Create jitter for x-positions
    jitter = 0.2 * randn(size(data, 1), 1); % Adjust 0.2 for spread
    
    % Plot individual points
    scatter(xPos + jitter, data(:, group), ...
        40, ... % Marker size
        boxColors(group,:), ... % Color matches box
        'filled', ...
        'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.5);
end
hold off;

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