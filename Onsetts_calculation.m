
clear; clc;clearvars
% list []condition labelss to reorganize (combine data)
step = {'max', 'normal'};
exo = {'exo', 'noexo'};
leg = {'left', 'right'};
%trialnum = [7, 6; 7, 7, 6, 0]; % Change as needed
%cumnum = cumsum(trialnum);

datafolder = 'mat'; % folder contatining trial data
%cd(datafolder)
choose_exo = ['Enter'];
for i = 1:length(exo)
    choose_exo = [choose_exo, int2str(i) 'for' exo(i) ','];
end

choose_step = ['Enter'];
for i = 1:length(step)
    choose_step = [choose_step, int2str(i) 'for' step(i) ','];
end

choose_leg = ['Enter'];
for i = 1:length(leg)
    choose_leg = [choose_leg, int2str(i) 'for' leg(i) ','];
end
exo_choose = inputdlg(strjoin(string(choose_exo)));
step_choose = inputdlg(strjoin(string(choose_step)));
leg_choose = inputdlg(strjoin(string(choose_leg)));

exo_choose = str2num(exo_choose{1});
step_choose = str2num(step_choose{1});
leg_choose = str2num(leg_choose{1});

files = dir;
filename = {};
indx = 1;
for i = 1:length(files)
    if files(i).bytes ~= 0 && contains(files(i).name,'.mat')
        filename{indx} = files(i).name;
        indx = indx+1;
    end
end
filenames = cell2mat(filename);
filenames = split(filenames,'.mat');

selected_files = {};
indx = 1;
for i = 1:length(filenames)
    indx_x = find(filenames{i,:} == '_');
    if contains(filenames{i,:},exo{exo_choose}) && contains(filenames{i,:},step{step_choose}) && contains(filenames{i,:},leg{leg_choose}) && (indx_x(2) - indx_x(1)-1 == length(exo{exo_choose}))
        selected_files{indx} = filenames{i,:};
        indx = indx+1;
    end
end
trialnum = length(selected_files);
%%

for tt=1:trialnum
    filename = [selected_files{tt} '.mat'];
    load(filename)
    if tt==1
        S = whos('-file',filename);
        varnames = [];
        for vv=1:length(S)
            varnames{vv} = S(vv).name;
        end
        % variable selection in order of preference
        vars = [9, 11, 6, 16, 18, 14, 22, ... % joint kinematics: leg (left to right, proximal to distal), then thorax
            10, 13, 7, 17, 19, 15, ... % joing moments
            2, 3, ... COM COP
            12, 20]; %
        %                 % clear S
    end
    trial = tt;
    Data(trial).Cond = [step{step_choose} '_' exo{exo_choose} '_' leg{leg_choose} int2str(tt)];
    for vv = 1:length(varnames)
        %eval(['Data(trial).' varnames{vv} '=  ' varnames{vv} '{1,1};']);
        if contains(varnames{vv},'Moment')
            eval(['Data(trial).' varnames{vv} '=  ' varnames{vv} '{1,1}(3:end-2,:);']);
        else
            eval(['Data(trial).' varnames{vv} '=  ' varnames{vv} '{1,1};']);
            %Data(trial).varnames{vv}(1:2,:) = [];
            %Data(trial).varnames{vv}(end-1:end,:) = [];
        end   
    end
end

Fs_mot = 100; % motion capture sampling frequency
Fs_load = 1000; % loadcell sampling frequency

%% converting signs, based on axis (in mocap) vs.  



%% plot


% We conclude that using first inflection point of hip flexion acc could
% work
%% Filtering

% Low-pass filtering
F_low = 4; % cutoff lower bound
[b,a] = butter(2, F_low/(Fs_mot/2), 'low');

for tt = 1:size(Data,2) % Adjust the range as needed
        if leg_choose == 1
            if ~isempty(Data(tt).LGRF)
                hippos = Data(tt).LGRF(:,2);
            else
                continue;
            end
            %disp(hippos)
        elseif leg_choose == 2
            if ~isempty(Data(tt).RGRF)
                hippos = Data(tt).RGRF(:,2);
            else
                continue;
            end
            %disp(hippos)
        end
            time = (1:length(hippos))/Fs_mot;
    
            hipvel = gradient(hippos,time); % Computing hip velocity
            hipvel = filtfilt(b,a,hipvel);
        
            hipacc = gradient(hipvel,time); % Computing hip acceleration
            hipacc = filtfilt(b,a,hipacc);
            
            Data(tt).hipvel = hipvel;
            Data(tt).hipacc = hipacc; % Saving the data of hip velocity & acceleration
end

close all; clc;

for tt=1:size(Data,2)
    if ~isempty(Data(tt).LGRF)
        data = Data(tt).LGRF(:,1);
        time = (1:length(data))/Fs_mot;
    
        subplot(3,3,tt)
        yyaxis left
        plot(time, data)
        yyaxis right
        plot(time, Data(tt).hipacc);
    end
end

figure;
for tt=1:size(Data,2)
    if ~isempty(Data(tt).RGRF)
        data = Data(tt).RGRF(:,1);
        time = (1:length(data))/Fs_mot;
    
        subplot(3,3,tt)
        yyaxis left
        plot(time, data)
        yyaxis right
        plot(time, Data(tt).hipacc);
    end
end

%% Onset detection
for tt=1:size(Data,2)
    hipacc = Data(tt).hipacc;
    m_win_hipacc = mean(abs(hipacc(1:50))); % Taking mean of the first 50 samples
    std_win_hipacc = std(abs(hipacc(1:50))); % Taking std of the first 50 samples
    for i = 1:length(hipacc)/50
        win_hipacc = hipacc((i-1)*50+1:i*50+1);
        onset_start = find(abs(win_hipacc) > abs(m_win_hipacc+3*std_win_hipacc)); % for a given window if th values are above the mean +- sd of the previous 50 samples
        if ~isempty(onset_start) && i >1
            onset = (i-1)*50+1+onset_start(1);
            Data(tt).onset = onset; % Onset determination
            break;
        else
            m_win_hipacc = mean(abs(hipacc((i-1)*50+1:i*50+1))); % If not compute mean & std for the previous 50 samples
            std_win_hipacc = std(abs(hipacc((i-1)*50+1:i*50+1)));
        end
    end

    figure % Plotting onset
    plot(hipacc)
    hold on
    plot([onset,onset],[min(hipacc),max(hipacc)],'b--')
end


% Initialize array to store onset values
onsets = zeros(1, size(Data, 2));

% Compute means of the first 10 data points for each trial
for tt = 1:size(Data, 2)
    if leg_choose == 1
        hippos = 1*(Data(tt).RHipAngle(:, 1));
        hippos1 = Data(tt).RHipAngle(:, 1);
        hippos = filtfilt(b, a, hippos1);
    elseif leg_choose == 2
        hippos = 1*(Data(tt).LHipAngle(:, 1));
        hippos1 = Data(tt).LHipAngle(:, 1);
        hippos = filtfilt(b, a, hippos1);
    end
    
    % Find the peak of the hip angle
    [~, peakIndex] = min(hippos);
    
    % Use the peak index as the onset value
    onsets(tt) = peakIndex;
    
    % Store the hip angle data
    Data(tt).HipAnkle1 = hippos;
    
    % Plot the trial
    figure;
    plot(hippos1);
    hold on;
    plot(peakIndex, hippos1(peakIndex), 'ro'); % Mark the peak on the plot
    xlabel('GRF');
    ylabel('Value');
    title(['Trial ', num2str(tt)]);
    hold off;
    pause(1); % Pause to display the plot before closing
    close; % Close the figure after displaying
    
end
% Subtract the mean of the first 10 data points from each column
hipposs = zeros(size(Data, 2), 571); % Initialize hipposs matrix
for tt = 1:size(Data, 2)
    means = mean(Data(tt).HipAnkle1(1:100));
    Data(tt).HipAnkle1 = Data(tt).HipAnkle1 - repmat(means, size(Data(tt).HipAnkle1, 1), 1);
    hipposs(tt,:) =   Data(tt).HipAnkle1(onsets(tt)-300:onsets(tt)+270); % based on onset take -75 and +500 points
    %hippos_stance_fro(tt,:)=hippos_stance_fro(onsets(tt)-170:onsets(tt)+500);
    %hippos_stance_sag(tt,:)=hippos_stance_sag(onsets(tt)-170:onsets(tt)+500);
end
figure;
plot(hipposs')

% Compute mean and standard deviation of hipposs
m_hipposs = mean(hipposs, 1);
std_hipposs = std(hipposs, 1);

% Plot ensemble average of the hip based on the onset


% Save the variables into a MAT file
%csave(fullfile(folder_path, 'max_noexo_right_sagittal.mat'), 'hipposs', 'm_hipposs', 'std_hipposs', 'onsets');


%folder_path = '/Users/mounikapasavula/Desktop/Gait study /Dryrun _experiment/Moments';

% Save the variables into a MAT file in the specified folder
%save(fullfile(folder_path, 'max_exo_right_1st.mat'), 'hipposs', 'm_hipposs', 'std_hipposs', 'onsets');

%%
%clear; clc
clearvars -except m_hipposs std_hipposs onsets

step = {'max', 'normal'};
exo = {'exo', 'noexo'};
leg = {'left', 'right'};
%trialnum = [7, 6; 7, 7, 6, 0]; % Change as needed
%cumnum = cumsum(trialnum);

datafolder = 'mat'; % folder contatining trial data
%cd(datafolder)
choose_exo = ['Enter'];
for i = 1:length(exo)
    choose_exo = [choose_exo, int2str(i) 'for' exo(i) ','];
end

choose_step = ['Enter'];
for i = 1:length(step)
    choose_step = [choose_step, int2str(i) 'for' step(i) ','];
end

choose_leg = ['Enter'];
for i = 1:length(leg)
    choose_leg = [choose_leg, int2str(i) 'for' leg(i) ','];
end
exo_choose = inputdlg(strjoin(string(choose_exo)));
step_choose = inputdlg(strjoin(string(choose_step)));
leg_choose = inputdlg(strjoin(string(choose_leg)));

exo_choose = str2num(exo_choose{1});
step_choose = str2num(step_choose{1});
leg_choose = str2num(leg_choose{1});

files = dir;
filename = {};
indx = 1;
for i = 1:length(files)
    if files(i).bytes ~= 0 && contains(files(i).name,'.mat')
        filename{indx} = files(i).name;
        indx = indx+1;
    end
end
filenames = cell2mat(filename);
filenames = split(filenames,'.mat');

selected_files = {};
indx = 1;
for i = 1:length(filenames)
    indx_x = find(filenames{i,:} == '_');
    if contains(filenames{i,:},exo{exo_choose}) && contains(filenames{i,:},step{step_choose}) && contains(filenames{i,:},leg{leg_choose}) && (indx_x(2) - indx_x(1)-1 == length(exo{exo_choose}))
        selected_files{indx} = filenames{i,:};
        indx = indx+1;
    end
end
trialnum = length(selected_files);
%%

for tt=1:trialnum
    filename = [selected_files{tt} '.mat'];
    load(filename)
    if tt==1
        S = whos('-file',filename);
        varnames = [];
        for vv=1:length(S)
            varnames{vv} = S(vv).name;
        end
        % variable selection in order of preference
        vars = [9, 11, 6, 16, 18, 14, 22, ... % joint kinematics: leg (left to right, proximal to distal), then thorax
            10, 13, 7, 17, 19, 15, ... % joing moments
            2, 3, ... COM COP
            12, 20]; %
        %                 % clear S
    end
    trial = tt;
    Data(trial).Cond = [step{step_choose} '_' exo{exo_choose} '_' leg{leg_choose} int2str(tt)];
    for vv = 1:length(varnames)
        %eval(['Data(trial).' varnames{vv} '=  ' varnames{vv} '{1,1};']);
        if contains(varnames{vv},'Moment')
            eval(['Data(trial).' varnames{vv} '=  ' varnames{vv} '{1,1}(3:end-2,:);']);
        else
            eval(['Data(trial).' varnames{vv} '=  ' varnames{vv} '{1,1};']);
            %Data(trial).varnames{vv}(1:2,:) = [];
            %Data(trial).varnames{vv}(end-1:end,:) = [];
        end   
    end
end

Fs_mot = 100; % motion capture sampling frequency
Fs_load = 1000; % loadcell sampling frequency

%% converting signs, based on axis (in mocap) vs.  



%% plot


% We conclude that using first inflection point of hip flexion acc could
% work
%% Filtering

% Low-pass filtering
F_low = 4; % cutoff lower bound
[b,a] = butter(2, F_low/(Fs_mot/2), 'low');

for tt = 1:size(Data,2) % Adjust the range as needed
        if leg_choose == 1
            if ~isempty(Data(tt).LGRF)
                hippos = Data(tt).LGRF(:,2);
            else
                continue;
            end
            %disp(hippos)
        elseif leg_choose == 2
            if ~isempty(Data(tt).RGRF)
                hippos = Data(tt).RGRF(:,2);
            else
                continue;
            end
            %disp(hippos)
        end
            time = (1:length(hippos))/Fs_mot;
    
            hipvel = gradient(hippos,time); % Computing hip velocity
            hipvel = filtfilt(b,a,hipvel);
        
            hipacc = gradient(hipvel,time); % Computing hip acceleration
            hipacc = filtfilt(b,a,hipacc);
            
            Data(tt).hipvel = hipvel;
            Data(tt).hipacc = hipacc; % Saving the data of hip velocity & acceleration
end

close all; clc;

for tt=1:size(Data,2)
    if ~isempty(Data(tt).LGRF)
        data = Data(tt).LGRF(:,1);
        time = (1:length(data))/Fs_mot;
    
        subplot(3,3,tt)
        yyaxis left
        plot(time, data)
        yyaxis right
        plot(time, Data(tt).hipacc);
    end
end

figure;
for tt=1:size(Data,2)
    if ~isempty(Data(tt).RGRF)
        data = Data(tt).RGRF(:,1);
        time = (1:length(data))/Fs_mot;
    
        subplot(3,3,tt)
        yyaxis left
        plot(time, data)
        yyaxis right
        plot(time, Data(tt).hipacc);
    end
end

%% Onset detection
for tt=1:size(Data,2)
    hipacc = Data(tt).hipacc;
    m_win_hipacc = mean(abs(hipacc(1:50))); % Taking mean of the first 50 samples
    std_win_hipacc = std(abs(hipacc(1:50))); % Taking std of the first 50 samples
    for i = 1:length(hipacc)/50
        win_hipacc = hipacc((i-1)*50+1:i*50+1);
        onset_start = find(abs(win_hipacc) > abs(m_win_hipacc+3*std_win_hipacc)); % for a given window if th values are above the mean +- sd of the previous 50 samples
        if ~isempty(onset_start) && i >1
            onset = (i-1)*50+1+onset_start(1);
            Data(tt).onset = onset; % Onset determination
            break;
        else
            m_win_hipacc = mean(abs(hipacc((i-1)*50+1:i*50+1))); % If not compute mean & std for the previous 50 samples
            std_win_hipacc = std(abs(hipacc((i-1)*50+1:i*50+1)));
        end
    end

    figure % Plotting onset
    plot(hipacc)
    hold on
    plot([onset,onset],[min(hipacc),max(hipacc)],'b--')
end


% Initialize array to store onset values
onsets1 = zeros(1, size(Data, 2));

% Compute means of the first 10 data points for each trial
for tt = 1:size(Data, 2)
    if leg_choose == 1
        hipos = 1*(Data(tt).RHipAngle(:, 1));
        hipos1 = Data(tt).RHipAngle(:, 1);
        hipos = filtfilt(b, a, hipos1);
    elseif leg_choose == 2
        hipos = 1*(Data(tt).LHipAngle(:, 1));
        hipos1 = Data(tt).LHipAngle(:, 1);
        hipos = filtfilt(b, a, hipos1);
    end
    
    % Find the peak of the hip angle
    [~, peakIndex1] = min(hipos);
    
    % Use the peak index as the onset value
    onsets1(tt) = peakIndex1;
    
    % Store the hip angle data
    Data(tt).HipAnkle11 = hipos;
    
    % Plot the trial
    figure;
    plot(hipos1);
    hold on;
    plot(peakIndex1, hipos1(peakIndex1), 'ro'); % Mark the peak on the plot
    xlabel('GRF');
    ylabel('Value');
    title(['Trial ', num2str(tt)]);
    hold off;
    pause(1); % Pause to display the plot before closing
    close; % Close the figure after displaying
    
end
% Subtract the mean of the first 10 data points from each column
hiposs = zeros(size(Data, 2), 571); % Initialize hipposs matrix
for tt = 1:size(Data, 2)
    means = mean(Data(tt).HipAnkle11(1:100));
    Data(tt).HipAnkle11 = Data(tt).HipAnkle11 - repmat(means, size(Data(tt).HipAnkle11, 1), 1);
    hiposs(tt,:) =   Data(tt).HipAnkle11(onsets1(tt)-300:onsets1(tt)+270); % based on onset take -75 and +500 points
    %hippos_stance_fro(tt,:)=hippos_stance_fro(onsets(tt)-170:onsets(tt)+500);
    %hippos_stance_sag(tt,:)=hippos_stance_sag(onsets(tt)-170:onsets(tt)+500);
end


% Compute mean and standard deviation of hipposs
m_hiposs = mean(hiposs, 1);
std_hiposs = std(hiposs, 1);

% Plot ensemble average of the hip based on the onset


% Save the variables into a MAT file
%csave(fullfile(folder_path, 'max_noexo_right_sagittal.mat'), 'hipposs', 'm_hipposs', 'std_hipposs', 'onsets');


%folder_path = '/Users/mounikapasavula/Desktop/Gait study /Dryrun _experiment/Moments';

% Save the variables into a MAT file in the specified folder
%save(fullfile(folder_path, 'max_exo_right_1st.mat'), 'hipposs', 'm_hipposs', 'std_hipposs', 'onsets');


%plot(m_hiposs, 'r-', 'LineWidth', 2);
%plot(m_hiposs - std_hiposs, 'r--');
%plot(m_hiposs + std_hiposs, 'r--');


% Plot the first plot with its standard deviation shaded in lighter red
% Plot the first plot with its standard deviation shaded in lighter red
figure;
%plot(m_hipposs, 'color',[83/255, 27/255, 147/255] ,'LineWidth', 3.5);
plot(m_hipposs, 'color',[148/255, 17/255, 0] ,'LineWidth', 3.5);
hold on;
x = 1:numel(m_hipposs);
%fill([x fliplr(x)], [m_hipposs - std_hipposs fliplr(m_hipposs + std_hipposs)], [83/255, 27/255, 147/255], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
fill([x fliplr(x)], [m_hipposs - std_hipposs fliplr(m_hipposs + std_hipposs)], [148/255, 17/255, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);


% Plot the second plot with its standard deviation shaded in lighter orange

%plot(m_hiposs, 'Color', [0, 84/255, 147/255], 'LineWidth', 3.5); % Orange color
plot(m_hiposs, 'Color', [255/255, 147/255, 0], 'LineWidth', 3.5);
hold on;
x = 1:numel(m_hiposs);
%fill([x fliplr(x)], [m_hiposs - std_hiposs fliplr(m_hiposs + std_hiposs)], [0, 84/255, 147/255], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
fill([x fliplr(x)], [m_hiposs - std_hiposs fliplr(m_hiposs + std_hiposs)], [255/255, 147/255, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
xlabel('Time (s)');
%ylabel('Adduction (-) / Abduction (+) Angles ( °)');
%ylabel('Assistance (-) / Resistance (+) Forces (N)');
ylabel('Adduction (-) / Abduction (+) Torques (N.m)');
%title('1');
xticks(100:100:700);
xticklabels({'1','2','3','4','5','6','7'});
yticks(-20:20:120);
ylim([-20, 120]);
hold off
% Save the variables into a MAT file
%save(fullfile(folder_path, 'max_noexo_right_sagittal.mat'), 'hipposs', 'm_hipposs', 'std_hipposs', 'onsets');


folder_path = '/Users/mounikapasavula/Documents/Gait study/All Subjects_Data';

% Save the variables into a MAT file in the specified folder
save(fullfile(folder_path, 'max_normal_noexo_right.mat'), 'm_hipposs', 'std_hipposs', 'onsets', 'm_hiposs', 'std_hiposs', 'onsets1');

