function tf_MVL_all = extract_and_analyze_eeg_time_resolved(start_time, end_time, filename, dataset_path, channel_name)
    % Function to perform time-resolved PAC analysis with sliding windows
    % Computes time-varying Phase-Amplitude Coupling (PAC) for each low frequency (2-12 Hz)
    % against the entire high frequency range (30-60 Hz) with a 100 ms step size.
    %
    % Inputs:
    % - start_time: Start time (in seconds) of the desired segment
    % - end_time: End time (in seconds) of the desired segment
    % - filename: Name of the EEG dataset (.set file)
    % - dataset_path: Folder path to the EEG dataset (.set file)
    % - channel_name: Channel to analyze (e.g., 'Fz')

    %% Necessary Paths
    warning off;
    addpath('functions');
    addpath('tfdnmfiles');
    addpath('EEG Data');
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

    %% Load Dataset
    EEG = pop_loadset(filename, dataset_path);
    EEG = eeg_checkset(EEG);
    Fs = EEG.srate;

    %% Channel Selection
    channel_idx = find(strcmp({EEG.chanlocs.labels}, channel_name));
    if isempty(channel_idx)
        error('Channel name not found in dataset. Check spelling or channel list.');
    end

    %% EEG Full Window Comparison
    EEG_Full_Segment = pop_select(EEG, 'time', [start_time end_time]);
    eeg_full_data = EEG_Full_Segment.data(channel_idx, :);
    
    %% Time Resolution
    window_length = 1; % 1 sec windows
    overlap = 0.9; % 90% overlap
    step_size = window_length * (1 - overlap); % Step size for sliding window
    total_windows = floor((end_time - start_time - window_length) / step_size) + 1;
    lowfreq = 2:12;
    highfreq = 30:60;
    tf_MVL_all = zeros(length(highfreq), length(lowfreq), total_windows);

    %% Compute TF-MVL Over Time
    for win = 1:total_windows
        win_start = start_time + (win - 1) * step_size;
        win_end = win_start + window_length;
        disp(win_start)
        disp(win_end)
        EEG_segment = pop_select(EEG, 'time', [win_start win_end]);
        eeg_data = EEG_segment.data(channel_idx, :);

        for i = 1:length(lowfreq)
            for j = 1:length(highfreq)
                tf_MVL_all(j, i, win) = tfMVL(eeg_data, highfreq(j), lowfreq(i), Fs);
            end
        end
    end
    disp(size(tf_MVL_all))
    
    % Plot Heatmaps in a Single Figure
figure;
num_plots = length(lowfreq) + 1;  % Total number of subplots

for i = 1:length(lowfreq)
    subplot(6, 2, i);  % Adjusts to fit in a 6x2 grid
    imagesc((1:total_windows) * 100, highfreq, squeeze(tf_MVL_all(:, i, :)));
    set(gca, 'YDir', 'normal');
    xlabel('Time (ms)');
    ylabel('High Frequency (Hz)');
    title(['Low Freq: ', num2str(lowfreq(i)), ' Hz']);
    colormap(jet);
    colorbar;
    axis tight;
end

% Add the EEG full data plot at the end
time_vector = linspace(0, (end_time - start_time) * 1000, length(eeg_full_data));
% time_vector = linspace(start_time * 1000, end_time * 1000, length(eeg_full_data));
subplot(6, 2, num_plots);  
plot(time_vector, eeg_full_data);
xlabel('Time (ms)');
ylabel('Amplitude');
title('EEG Full Data');

sgtitle('Time-Resolved PAC Heatmaps');
