% Clear workspace and command window
clear;
clc;

% Define muscle list with right and left indices for 16 channels
params.muslist = [
    3, 1; % Gmax right, left (2148 Hz)
    4, 2; % Gmed right, left (2148 Hz)
    6, 5; % Adl right, left (2222 Hz)
    8, 7; % RF right, left (2222 Hz)
    10, 9; % VL right, left (2222 Hz)
    14, 13; % BF right, left (2222 Hz)
    12, 11; % TA right, left (2148 Hz)
    16, 15; % MG right, left (2148 Hz)
];

% Define conditions
params.exo = {'BA', 'noexo', 'noexo', 'exo', 'exo'};
params.leg = {'N/A', 'left', 'left', 'left', 'left'};
params.step = {'N/A', 'normal', 'max', 'normal', 'max'};

% Define trial numbers for each condition
params.trialnum = {
    1:3, ... % Baseline trials
    4:10, ... % WO R Norm
    11:17, ... % W0 R Long
    18:24, ... % W R Norm
    25:31, ... % W R Long
};

% Initialize data storage structure
S01 = struct();
params.analysis_path = '/Users/mounikapasavula/Documents/Gait study/S02/';
params.EMG_path = '/Users/mounikapasavula/Documents/Gait study/S02/EMG_Mat/';
params.Torq_path = '/Users/mounikapasavula/Documents/Gait study/S02/QTM_mat/';

% Define target sample rate for data processing
params.EMG_fs = 2148; % 2148 Hz for channels 1-4 and 13-16
params.Torq_fs = 100;

% Initialize the trial index
trial_index = 1;

cd(params.EMG_path);
% Process each trial based on the conditions
for condition_index = 1:numel(params.trialnum)
    % Get the trials for the current condition
    trials = params.trialnum{condition_index};
    
    % Iterate through each trial
    for tt = 1:length(trials)
        % Construct the filename
        filename = sprintf('Run_number_32_Plot_and_Store_Rep_%d.%d.mat', trials(tt), trials(tt) - 1);
        
        % Load the data from the file
        load(filename);
        
        
        % Process each muscle channel for the current trial
        for i = 1:size(params.muslist, 1)
            channel_r = params.muslist(i, 1);
            channel_l = params.muslist(i, 2);
            
            % Get the sampling rates for each channel
            fs_r = round(Fs(channel_r)); % Sampling rate for right channel
            fs_l = round(Fs(channel_l)); % Sampling rate for left channel

            time_r = round(Time(channel_r,:),3); % Time for right channel
            time_l = round(Time(channel_l,:),3); % Time for left channel
            
            indx_r = find(time_r == 12);
            indx_l = find(time_l == 12);

            if length(indx_r) < 1
                indx_r = indx_r(1);
            end
            if length(indx_l) < 1
                indx_l = indx_l(1);
            end

            emg_r{i, :} = Data(channel_r, 1:indx_r);
            emg_l{i, :} = Data(channel_l, 1:indx_l);
            trial_time{i,:} = 0:12 / (length(emg_l{i,:}) - 1):12;
           
        end
        S01(trial_index).Exo = params.exo{condition_index};
        S01(trial_index).Leg = params.leg{condition_index};
        S01(trial_index).Step = params.step{condition_index};


        S01(trial_index).Right = emg_r;
        S01(trial_index).Left = emg_l;
        S01(trial_index).Time = trial_time;
        %save('S01.mat', 'S01');
        
        % Increment the trial index for the next trial
        trial_index = trial_index + 1;
        
        % Clear variables to free memory and prevent data leakage
        clear Data Time Fs;
    end
end
clearvars -except S01 params

cd(params.analysis_path)
%% Importing filenames for Torque

cd(params.Torq_path)
files = dir;
filenames = {};
indx = 0;
for i = 1:length(files)
    if files(i).bytes ~= 0 && files(i).isdir == 0 && (contains(files(i).name, '.mat') == 1)
        indx = indx + 1;
        filenames{indx} = files(i).name;
    end

end

%% combining torques
for i = 1:length(filenames)
    file = filenames{i};
    for j = 1:length(S01)
        if (contains(file,S01(j).Exo)==1) && (contains(file,S01(j).Leg)==1) && (contains(file,S01(j).Step)==1)
            load(file);
            if S01(j).Leg == 'left'
                S01(j).Hip_torque = RHipMoment{1}(3:end-2,:);
            elseif S01(j).Leg == 'right'
                S01(j).Hip_torque = LHipMoment{1}(3:end-2,:);
            end
        end
    end
end

%save('S01.mat', 'S01');
clearvars -except S01 params

cd(params.analysis_path)