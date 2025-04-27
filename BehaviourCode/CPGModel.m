% Enhanced Physiological CPG Model with Variable Population Dynamics
% References: 
% 1. Leaky Integrate-and-Fire (LIF) neuron model (Dayan & Abbott, 2001)
% 2. CPG-inspired reciprocal inhibition (Marder & Bucher, 2001)
% 3. Trappenberg (2010) for population modeling
% 4. Harris-Warrick (2011) for neuronal heterogeneity in CPGs
% 5. Grillner & El Manira (2020) for stochastic properties of CPGs

clc;

%% Parameters
n_neurons = 50;      % Total neurons (half flexors, half extensors)
dt = 0.1;            % Time step (ms)
t_total = 10000;      % Simulation time (ms)
V_rest = -65;        % Resting potential (mV)
V_reset = -70;       % Reset potential (mV)
R_m = 10;            % Membrane resistance (MΩ)
I_app = 2.2;         % Base input current (nA)
noise_amp = 0.3;     % Input noise amplitude

% Network parameters
g_syn = 0.4;         % Synaptic conductance
E_syn = -80;         % Inhibitory reversal potential
tau_syn = 5;         % Synaptic decay (ms)

% Pulsed Input Parameters (insert near other parameters)
pulse_interval = 1000;    % Input period (ms)
pulse_duration = 500;      % Active phase duration (ms)
I_app_pulse = 3.2;       % Increased amplitude for pulsed operation

%%% Physiological Enhancements - 1: Heterogeneous Neuron Properties
% Distribute membrane time constants
tau_m_mean = 10;     % Mean membrane time constant (ms)
tau_m_var = 1;       % Variance
tau_m = tau_m_mean + tau_m_var * randn(n_neurons, 1);
tau_m = max(tau_m, 5);  % Ensure values don't go too low
tau_pulse_integration = 8;  % Membrane integration time constant (ms)
% Distribute thresholds
V_thresh_mean = -50; % Mean threshold (mV)
V_thresh_var = 3;    % Variance (mV)
V_thresh = V_thresh_mean + V_thresh_var * randn(n_neurons, 1);

%%% Physiological Enhancements - 2: Stochastic Synaptic Transmission
syn_release_prob = 0.90;  % Probability of synaptic release

%%% Physiological Enhancements - 3: Heterogeneous Connectivity
% Create probabilistic connection matrix instead of all-to-all
conn_prob = 0.7;  % 70% connection probability
conn_matrix = rand(n_neurons, n_neurons) < conn_prob;

%%% Physiological Enhancements - 4: Refractory Period
ref_period = 3;      % Refractory period (ms)
ref_timer = zeros(n_neurons, 1);  % Refractory timers for each neuron

%% Initialize
t = 0:dt:t_total;
flexor_neurons = 1:n_neurons/2;
extensor_neurons = (n_neurons/2+1):n_neurons;

% Membrane potentials and spikes (neurons x time)
V = V_rest * ones(n_neurons, length(t));
spikes = zeros(n_neurons, length(t));

% Add some initial variability to membrane potentials
V(:,1) = V_rest + 1*randn(n_neurons, 1);

% Synaptic currents
I_syn = zeros(n_neurons, 1);

%% Simulation
for i = 2:length(t)
    % Physiological Enhancements - 5: Continuous Background Noise
    background_noise = 0.005 * randn(n_neurons, 1);  % 0.15 nA SD continuous noise
    % Modified Input Section (replace existing I_noise calculation)
    % In simulation loop (before neuron updates):
    current_time = t(i);
    if mod(current_time, pulse_interval) < pulse_duration
        I_app_current = I_app_pulse;
    else
        I_app_current = I_app;
    end
    % Revised noise/input calculation
    I_noise = I_app_current + noise_amp * randn(n_neurons, 1) + background_noise;

    % Update all neurons
    for n = 1:n_neurons
        % Skip if neuron is in refractory period
        if ref_timer(n) > 0
            ref_timer(n) = ref_timer(n) - dt;
            continue;
        end
        
        % Modify voltage equation:
        dV = (-(V(n,i-1) - V_rest) + R_m*(I_noise(n) - I_syn(n))) / ...
            (tau_m(n) + tau_pulse_integration*(I_app_current>0)) * dt;
        V(n,i) = V(n,i-1) + dV;

        % Physiological Enhancements - 6: Spontaneous Activity
        if rand < 0.0005  % Small probability of spontaneous spike
            V(n,i) = V_thresh(n) + 1;  % Force spike
        end
        
        % Spike detection
        if V(n,i) >= V_thresh(n)
            spikes(n,i) = 1;
            V(n,i) = V_reset;
            ref_timer(n) = ref_period;  % Set refractory period
            
            % Trigger synaptic inhibition to opposing population
            if ismember(n, flexor_neurons)
                targets = extensor_neurons;
                % Get connected targets based on connectivity matrix
                connected_targets = targets(conn_matrix(n, targets) == 1);
                % Apply stochastic synaptic transmission
                active_targets = connected_targets(rand(size(connected_targets)) < syn_release_prob);
                
                if ~isempty(active_targets)
                    % Add variability to synaptic strength (±20%)
                    g_effective = 0.25 * g_syn * (0.8 + 0.4*rand());
                    I_syn(active_targets) = g_effective * (V(active_targets,i-1) - E_syn);
                end
            else
                targets = flexor_neurons;
                % Get connected targets based on connectivity matrix
                connected_targets = targets(conn_matrix(n, targets) == 1);
                % Apply stochastic synaptic transmission
                active_targets = connected_targets(rand(size(connected_targets)) < syn_release_prob);
                
                if ~isempty(active_targets)
                    % Add variability to synaptic strength (±20%)
                    g_effective = 0.25 * g_syn * (0.8 + 0.4*rand());
                    I_syn(active_targets) = g_effective * (V(active_targets,i-1) - E_syn);
                end
            end
        end
    end
    
    % Decay synaptic currents
    I_syn = I_syn * exp(-dt/tau_syn);
end

%% Rastergram Plot
figure;
hold on;
title('CPG Network Rastergram');

% Flexor spikes (blue)
[f_neurons, f_times] = find(spikes(flexor_neurons,:) == 1);
scatter(t(f_times), f_neurons, 10, 'b', 'filled');

% Extensor spikes (red)
[e_neurons, e_times] = find(spikes(extensor_neurons,:) == 1);
scatter(t(e_times), e_neurons + n_neurons/2, 10, 'r', 'filled');

xlabel('Time (ms)'); 
ylabel('Neuron ID');
ylim([0 n_neurons+1]);
legend('Flexors', 'Extensors');

%% Population Activity Plot
figure;
subplot(2,1,1);
plot(t, mean(V(flexor_neurons,:)), 'b');
hold on;
plot(t, mean(V(extensor_neurons,:)), 'r');
title('Population Membrane Potential');
legend('Flexors', 'Extensors');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');

subplot(2,1,2);
flexor_rate = movmean(sum(spikes(flexor_neurons,:)), 50)/length(flexor_neurons)*1000/dt;
extensor_rate = movmean(sum(spikes(extensor_neurons,:)), 50)/length(extensor_neurons)*1000/dt;
plot(t, flexor_rate, 'b', t, extensor_rate, 'r');
title('Population Firing Rate (Hz)');
xlabel('Time (ms)');
ylabel('Firing Rate (Hz)');

%% Calculate phase relationships

% Phase coherence analysis
figure;

% Extract burst times from population activity
flexor_pop = sum(spikes(flexor_neurons,:), 1);
extensor_pop = sum(spikes(extensor_neurons,:), 1);

% Find peaks in population activity (bursts)
[~, flexor_peaks] = findpeaks(flexor_rate, 'MinPeakHeight', 50, 'MinPeakDistance', 20);
[~, extensor_peaks] = findpeaks(extensor_rate, 'MinPeakHeight', 50, 'MinPeakDistance', 20);

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
polarhistogram(phase_diff, 20, 'Normalization', 'probability');
title('Phase Relationship Between Flexor and Extensor Populations');

% Calculate phase coherence metric (R)
R = abs(mean(exp(1i*phase_diff)));
phase_coherence = R * 100; % Convert to percentage
disp(['Phase coherence: ', num2str(phase_coherence), '%']);
%%
% Calculate locomotion metrics from CPG output
figure;

% 1. Stride length (proportional to cycle duration)
mean_cycle_length = mean(cycle_lengths);
std_cycle_length = std(cycle_lengths);

% 2. Normalize cycle lengths to a reference value (for mice ~200ms)
reference_stride = 200; % ms
stride_lengths = cycle_lengths * (reference_stride/mean_cycle_length);

% 3. Calculate stride variability metrics
CV_stride = std_cycle_length / mean_cycle_length * 100; % Coefficient of variation (%)

% 4. Plot stride length distribution
subplot(2,1,1);
histogram(stride_lengths, 10);
title('Stride Length Distribution');
xlabel('Stride Length (normalized)');
ylabel('Count');

% 5. Plot stride-to-stride variability
subplot(2,1,2);
plot(1:length(stride_lengths), stride_lengths, 'o-');
title(['Stride-to-Stride Variability (CV = ', num2str(CV_stride), '%)']);
xlabel('Stride Number');
ylabel('Stride Length');

% Display metrics
disp(['Mean stride length: ', num2str(mean(stride_lengths)), ' (normalized)']);
disp(['Stride length CV: ', num2str(CV_stride), '%']);
%%
% Calculate swing time metrics from CPG output
figure;

% Extract burst times from population activity
flexor_pop = sum(spikes(flexor_neurons,:), 1);
extensor_pop = sum(spikes(extensor_neurons,:), 1);

% Smooth the population activity
smoothed_flexor = movmean(flexor_pop, 20);
smoothed_extensor = movmean(extensor_pop, 20);

% Set activity thresholds
flexor_threshold = max(smoothed_flexor) * 0.25;
extensor_threshold = max(smoothed_flexor) * 0.25;

% Identify swing phases (when flexor activity > threshold)
swing_active = smoothed_flexor > flexor_threshold;
swing_phase_starts = find(diff([0 swing_active]) == 1);
swing_phase_ends = find(diff([swing_active 0]) == -1);

% Calculate swing time durations
swing_times = (swing_phase_ends - swing_phase_starts) * dt*1000;

% Plot swing time distribution and variability
subplot(2,1,1);
histogram(swing_times, 10);
title('Swing Time Distribution');
xlabel('Swing Time (ms)');
ylabel('Count');

subplot(2,1,2);
plot(1:length(swing_times), swing_times, 'o-b');
hold on;
plot([1 length(swing_times)], [mean(swing_times) mean(swing_times)], '--k');
title(['Swing Time (Mean: ', num2str(mean(swing_times)), ' ms, CV: ', num2str(std(swing_times)/mean(swing_times)*100), '%)']);
xlabel('Cycle Number');
ylabel('Swing Time (ms)');

% Display metrics
disp(['Mean swing time: ', num2str(mean(swing_times)), ' ms']);
disp(['Swing time variability (SD): ', num2str(std(swing_times)), ' ms']);
disp(['Coefficient of variation: ', num2str(std(swing_times)/mean(swing_times)*100), '%']);

% Calculate swing/stance ratio (important gait parameter)
stance_times = [];
for i = 1:length(swing_phase_ends)-1
    stance_times(i) = (swing_phase_starts(i+1) - swing_phase_ends(i)) * dt;
end

if ~isempty(stance_times)
    swing_stance_ratio = mean(swing_times) / mean(stance_times);
    disp(['Swing/stance ratio: ', num2str(swing_stance_ratio)]);
end
