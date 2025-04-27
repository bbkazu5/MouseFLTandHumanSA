%% Create publication-quality figure with multiple panels
figure('Position', [50 50 1000 800], 'Color', 'w');  % Large white figure
tl = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Panel A: Population Membrane Potential (moved from Panel D)
nexttile([1 2]);
plot(t, mean(V(flexor_neurons,:)), 'b', 'LineWidth', 1.5);
hold on;
plot(t, mean(V(extensor_neurons,:)), 'r', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('A: Population Membrane Potential', 'FontWeight', 'normal');
legend('Flexors', 'Extensors', 'Location', 'northeastoutside');
box on;
xlim([100 650]);
%% Panel B: Rastergram of neural activity
nexttile(4, [1 2]);
hold on;
[f_neurons, f_times] = find(spikes(flexor_neurons,:) == 1);
scatter(t(f_times), f_neurons, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.7);
[e_neurons, e_times] = find(spikes(extensor_neurons,:) == 1);
scatter(t(e_times), e_neurons + n_neurons/2, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.7);
xlim([100 650]);  % Zoom to a representative section
xlabel('Time (ms)'); 
ylabel('Neuron ID');
ylim([0 n_neurons+1]);
title('B: CPG Network Activity', 'FontWeight', 'normal');
legend('Flexors', 'Extensors', 'Location', 'northeastoutside');
box on;

%% Panel C: Smoothed Firing Rates (for muscle contraction)
nexttile(7, [1 2]);
% Calculate population firing rates with stronger smoothing for muscle contraction
window_size = 20;  % Larger window for smoother muscle-like activation
flexor_smooth = movmean(sum(spikes(flexor_neurons,:), 1), window_size)/length(flexor_neurons)*1000/dt;
extensor_smooth = movmean(sum(spikes(extensor_neurons,:), 1), window_size)/length(extensor_neurons)*1000/dt;

% Normalize to 0-1 range for muscle activation interpretation
flexor_muscle = flexor_smooth/max(flexor_smooth);
extensor_muscle = extensor_smooth/max(flexor_smooth); % Use same scale as flexor

plot(t, flexor_muscle, 'b', 'LineWidth', 2);
hold on;
plot(t, extensor_muscle, 'r', 'LineWidth', 2);
xlim([100 650]);  % Match the same time window as panels above
xlabel('Time (ms)');
ylabel('Normalized Muscle Activation');
title('C: Muscle Activation Patterns', 'FontWeight', 'normal');
box on;

%% Panel E: Phase Relationship Polar Plot
nexttile(3, [1 1]);
% Extract burst times from population activity
flexor_pop = sum(spikes(flexor_neurons,:), 1);
extensor_pop = sum(spikes(extensor_neurons,:), 1);

% Find peaks in population activity (bursts)
[~, flexor_peaks] = findpeaks(flexor_smooth, 'MinPeakHeight', max(flexor_smooth)*0.2, 'MinPeakDistance', 20);
[~, extensor_peaks] = findpeaks(extensor_smooth, 'MinPeakHeight', max(extensor_smooth)*0.2, 'MinPeakDistance', 20);

% Calculate phase differences (in radians)
phase_diff = zeros(length(flexor_peaks), 1);
cycle_lengths = zeros(length(flexor_peaks)-1, 1);

for i = 1:length(flexor_peaks)-1
    % Cycle length is time between consecutive flexor bursts
    cycle_length = (flexor_peaks(i+1) - flexor_peaks(i)) * dt;
    cycle_lengths(i) = cycle_length;
    
    % Find nearest extensor peak within this cycle
    ext_in_cycle = extensor_peaks(extensor_peaks > flexor_peaks(i) & extensor_peaks < flexor_peaks(i+1));
    
    if ~isempty(ext_in_cycle)
        % Calculate phase as fraction of cycle (0-2π)
        phase = 2*pi * (ext_in_cycle(1) - flexor_peaks(i)) * dt / cycle_length;
        phase_diff(i) = phase;
    end
end

% Remove zeros (cycles without extensor peaks)
phase_diff = phase_diff(phase_diff > 0);

% Create polar plot of phase differences
polarhistogram(phase_diff, 20, 'Normalization', 'probability', 'FaceColor', [0.6 0.6 0.9]);
title('E: Phase Relationship', 'FontWeight', 'normal');

% Calculate phase coherence metric (R)
R = abs(mean(exp(1i*phase_diff)));
phase_coherence = R * 100; % Convert to percentage
text(0, 0, ['Coherence: ', num2str(round(phase_coherence)), '%'], 'HorizontalAlignment', 'center');
%% Panel F: Stride Length Analysis
nexttile(6, [1 1]);
% 1. Stride length (proportional to cycle duration)
mean_cycle_length = mean(cycle_lengths);
std_cycle_length = std(cycle_lengths);

% 2. Normalize cycle lengths to a reference value (for mice ~200ms)
reference_stride = 100; % ms
stride_lengths = cycle_lengths * (reference_stride/mean_cycle_length);

% 3. Calculate stride variability metrics
CV_stride = std_cycle_length / mean_cycle_length * 100; % Coefficient of variation (%)

histogram(stride_lengths, 10, 'FaceColor', [0.4 0.7 0.9]);
title('F: Stride Length Distribution', 'FontWeight', 'normal');
xlabel('Stride Length (norm.)');
ylabel('Count');
text(min(stride_lengths)+range(stride_lengths)*0.1, max(histcounts(stride_lengths,10))*0.8,...
    ['CV = ', num2str(round(CV_stride,1)), '%'], 'FontSize', 9);
box on;
axis square
%% Panel G: Swing Time Distribution (fixed to ensure it's included)
nexttile(9, [1 1]);
% Process flexor activity to identify swing phases
flexor_threshold = max(flexor_smooth) * 0.25;

% Identify swing phases (when flexor activity > threshold)
swing_active = flexor_smooth > flexor_threshold;
swing_phase_starts = find(diff([0 swing_active]) == 1);
swing_phase_ends = find(diff([swing_active 0]) == -1);

% Calculate swing time durations
swing_times = (swing_phase_ends - swing_phase_starts) * dt;

% Plot swing time distribution
histogram(swing_times, 10, 'FaceColor', [0.9 0.6 0.4]);
title('G: Swing Time Distribution', 'FontWeight', 'normal');
xlabel('Swing Time (ms)');
ylabel('Count');
swing_CV = std(swing_times)/mean(swing_times)*100;
text(min(swing_times)+range(swing_times)*0.1, max(histcounts(swing_times,10))*0.8,...
    ['CV = ', num2str(round(swing_CV,1)), '%'], 'FontSize', 9);
box on;
axis square
%% Panel H: Swing vs Stride Relationship
nexttile(12, [1 1]);
% Calculate stance times
stance_times = [];
for i = 1:length(swing_phase_ends)-1
    if swing_phase_starts(i+1) > swing_phase_ends(i)  % ensure valid sequence
        stance_times(end+1) = (swing_phase_starts(i+1) - swing_phase_ends(i)) * dt;
    end
end

% Match stance and swing times to create pairs
min_length = min(length(stance_times), length(swing_times));
if min_length > 0
    % Calculate swing/stance ratio for each cycle
    ratio_values = swing_times(1:min_length) ./ stance_times(1:min_length);
    % Stride = swing + stance
    stride_from_phases = swing_times(1:min_length) + stance_times(1:min_length);
    
    % Plot relationship between swing time and stride length
    scatter(swing_times(1:min_length), stride_from_phases, 40, [0.4 0.5 0.8], 'filled', 'MarkerFaceAlpha', 0.7);
    xlabel('Swing Time (ms)');
    ylabel('Stride Length (ms)');
    title('H: Swing Time vs Stride Length', 'FontWeight', 'normal');
    
    % Add regression line
    hold on;
    p = polyfit(swing_times(1:min_length), stride_from_phases, 1);
    x_range = linspace(min(swing_times), max(swing_times), 100);
    y_fit = polyval(p, x_range);
    plot(x_range, y_fit, 'k--');
    R = corrcoef(swing_times(1:min_length), stride_from_phases);
    text(min(swing_times)+range(swing_times)*0.1, max(stride_from_phases)*0.9,...
        ['r = ', num2str(round(R(1,2),2))], 'FontSize', 9);
    box on;
end
axis square
%% Panel I: Swing/Stance Ratio Over Time
nexttile(10, [1 2]);
if min_length > 0
    % Plot swing/stance ratio over time (for cycles)
    cycle_times = swing_phase_starts(1:min_length) * dt; % Approximate cycle start times
    plot(cycle_times, ratio_values, 'o-', 'LineWidth', 1.5, 'Color', [0.8 0.4 0.2], 'MarkerSize', 5);
    xlabel('Time (ms)');
    ylabel('Swing/Stance Ratio');
    title('I: Swing/Stance Ratio', 'FontWeight', 'normal');
    % Add reference line for balanced gait
    hold on;
    yline(1, 'k--', 'Balanced Gait');
    box on;
end
xlim([100 650])
%% Overall title and formatting
title(tl, 'Central Pattern Generator Network and Derived Locomotion Metrics', 'FontSize', 14);

