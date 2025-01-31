
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

%close all; clc;
%{
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
%}
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

    %figure % Plotting onset
    %plot(hipacc)
    %hold on
    %plot([onset,onset],[min(hipacc),max(hipacc)],'b--')
end


% Initialize array to store onset values

load('a_noecu.mat', 'onsets');

% Compute means of the first 10 data points for each trial
for tt = 1:size(Data, 2)
    if leg_choose == 1
        hippos = 1 * (Data(tt).RHipMoment(:, 1));
        hippos1 = Data(tt).RGRF(:, 1);
        hippos = filtfilt(b, a, hippos);
    elseif leg_choose == 2
        hippos = -1* (Data(tt).LHipMoment(:, 1));
        %hippos=(5.0008*(hippos_) + 0.0028211)*4.44822161533;
        hippos1 = Data(tt).LGRF(:, 1);
        hippos = filtfilt(b, a, hippos);
    end
    
    % Use the pre-defined onset from the onsets variable
    onset_value = onsets(tt);
    %{
    % Plot the trial
    figure;
    plot(hippos1);
    xlabel('GRF');
    ylabel('Value');
    title(['Trial ', num2str(tt)]);
    
    % Mark the onset on the plot
    hold on;
    plot(onset_value, hippos1(onset_value), 'ro'); % Mark the onset with a red circle
    hold off;
    %}
    
    % Pause to visually confirm the onset (optional)
    %pause(1); % Adjust the pause duration as needed
    %close; % Close the figure after displaying
    
    Data(tt).HipAnkle1 = hippos;
end
% Subtract the mean of the first 10 data points from each column
hipposs = zeros(size(Data, 2), 571); % Initialize hipposs matrix
for tt = 1:size(Data, 2)
    means = mean(Data(tt).HipAnkle1(1:100));
    Data(tt).HipAnkle1 = Data(tt).HipAnkle1 - repmat(means, size(Data(tt).HipAnkle1, 1), 1);
    hipposs(tt,:) =   Data(tt).HipAnkle1(onsets(tt)-200:onsets(tt)+370); % based on onset take -75 and +500 points
    %hippos_stance_fro(tt,:)=hippos_stance_fro(onsets(tt)-170:onsets(tt)+500);
    %hippos_stance_sag(tt,:)=hippos_stance_sag(onsets(tt)-170:onsets(tt)+500);
end
%figure;
%plot(hipposs')

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

%close all; clc;
%{
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
%}
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

    %figure % Plotting onset
    %plot(hipacc)
    %hold on
    %plot([onset,onset],[min(hipacc),max(hipacc)],'b--')
end


% Initialize array to store onset values
load('a_noecu.mat', 'onsets1');

% Compute means of the first 10 data points for each trial
for tt = 1:size(Data, 2)
    if leg_choose == 1
        hipos = 1*(Data(tt).RHipMoment(:, 2));
        hipos1 = Data(tt).RGRF(:, 1);
        %hipos = (4.9984*(hipos) + 0.0062595)*4.4482216153;
        %hippos_stance_fro=Data(tt).RHipAngle(:, 2);
        %hippos_stance_sag=filtfilt(b, a, hippos_stance_sag);
        %hippos_stance_fro=filtfilt(b, a, hippos_stance_fro);
        hipos = filtfilt(b, a, hipos);
    elseif leg_choose == 2
        hipos = -1*(Data(tt).LHipMoment(:, 1));
        %hipos= (5.0008*(hipos_) + 0.0028211)*4.44822161533;
        hipos1 = Data(tt).LGRF(:, 1);
        %hipos = (5.0008*(hipos) + 0.0028211)*4.4482216153;
        %hippos_stance_sag=Data(tt).LHipAnkle(:, 1);
        %hippos_stance_fro=Data(tt).LHipAnkle(:, 2);
        %%hippos_stance_fro=filtfilt(b, a, hippos_stance_fro);
        hipos = filtfilt(b, a, hipos);
    end
    % Use the pre-defined onset from the onsets variable
    onset_value1 = onsets1(tt);

    %{
    % Plot the trial
    figure;
    plot(hipos1);
    xlabel('GRF');
    ylabel('Value');
    title(['Trial ', num2str(tt)]);
    
    % Mark the onset on the plot
    hold on;
    plot(onset_value1, hipos1(onset_value1), 'ro'); % Mark the onset with a red circle
    hold off;
    
    % Pause to visually confirm the onset (optional)
    pause(1); % Adjust the pause duration as needed
    close; % Close the figure after displaying
    %}
    
    Data(tt).HipAnkle11 = hipos;
end
% Subtract the mean of the first 10 data points from each column
hiposs = zeros(size(Data, 2), 571); % Initialize hipposs matrix
for tt = 1:size(Data, 2)
    means = mean(Data(tt).HipAnkle11(1:100));
    Data(tt).HipAnkle11 = Data(tt).HipAnkle11 - repmat(means, size(Data(tt).HipAnkle11, 1), 1);
    hiposs(tt,:) =   Data(tt).HipAnkle11(onsets1(tt)-200:onsets1(tt)+370); % based on onset take -75 and +500 points
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

subplot(2,2,1)
plot(m_hipposs, 'color',[83/255, 27/255, 147/255] ,'LineWidth', 3.5);
%plot(m_hipposs, 'color',[148/255, 17/255, 0] ,'LineWidth', 3.5);
hold on;
x = 1:numel(m_hipposs);
fill([x fliplr(x)], [m_hipposs - std_hipposs fliplr(m_hipposs + std_hipposs)], [83/255, 27/255, 147/255], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
%fill([x fliplr(x)], [m_hipposs - std_hipposs fliplr(m_hipposs + std_hipposs)], [148/255, 17/255, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);


% Plot the second plot with its standard deviation shaded in lighter orange

plot(m_hiposs, 'Color', [0, 84/255, 147/255], 'LineWidth', 3.5); % Orange color
%plot(m_hiposs, 'Color', [255/255, 147/255, 0], 'LineWidth', 3.5);
hold on;
x = 1:numel(m_hiposs);
fill([x fliplr(x)], [m_hiposs - std_hiposs fliplr(m_hiposs + std_hiposs)], [0, 84/255, 147/255], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
%fill([x fliplr(x)], [m_hiposs - std_hiposs fliplr(m_hiposs + std_hiposs)], [255/255, 147/255, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
xlabel('Time (s)');
%ylabel('Adduction (-) / Abduction (+) Angles ( °)');
%ylabel('Assistance (-) / Resistance (+) Forces (N)');
%ylabel('Adduction (-) / Abduction (+) Torques (N.m)');
%title('1');
xticks(100:100:700);
xticklabels({'1','2','3','4','5','6','7'});
yticks(-80:20:40);
ylim([-80, 40]);
fontsize(16,"points");
hold off
% Save the variables into a MAT file
%save(fullfile(f0]older_path, 'max_noexo_right_sagittal.mat'), 'hipposs', 'm_hipposs', 'std_hipposs', 'onsets');


%folder_path = '/Users/mounikapasavula/Documents/Gait study /S01/StanceAngles';

% Save the variables into a MAT file in the specified folder
%save(fullfile(folder_path, 'max_normal_noexo_right_Offset.mat'), 'm_hipposs', 'std_hipposs', 'onsets', 'm_hiposs', 'std_hiposs', 'onsets1');

%% Repeating for 2nd time 
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

%{
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

%}
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

    %figure % Plotting onset
    %plot(hipacc)
    %hold on
    %plot([onset,onset],[min(hipacc),max(hipacc)],'b--')
end


% Initialize array to store onset values

load('a_ecu.mat', 'onsets');

% Compute means of the first 10 data points for each trial
% Compute means of the first 10 data points for each trial
for tt = 1:size(Data, 2)
    if leg_choose == 1
        hippos = 1*(Data(tt).RHipMoment(:, 1));
        hippos1 = Data(tt).RGRF(:, 1);
        %hippos = (4.9984*(hippos) + 0.0062595)*4.4482216153;
        %hippos_stance_fro=Data(tt).RHipAngle(:, 2);
        %hippos_stance_sag=filtfilt(b, a, hippos_stance_sag);
        %hippos_stance_fro=filtfilt(b, a, hippos_stance_fro);
        hippos = filtfilt(b, a, hippos);
    elseif leg_choose == 2
        loadcell11 = 1*(Data(tt).LLoadCell);
        loadcell = (5.0008*(loadcell11) + 0.0028211)*4.44822161533;
        
        RB1_1 = [(Data(tt).LB1(:, 1)), (Data(tt).LB1(:, 2)), (Data(tt).LB1(:, 3))];
        RB2_1 = [(Data(tt).LB2(:, 1)), (Data(tt).LB2(:, 2)), (Data(tt).LB2(:, 3))];
        RFT_1= [Data(tt).LFT(:, 1), Data(tt).LFT(:, 2), Data(tt).LFT(:, 3)];
        
        % Calculate the displacement vector from RB1 to RB2
        displacement__vector = RB1_1 - RB2_1;
        displacement_vector1 = RFT_1;%- RB1_1;
        forceV= displacement__vector .* loadcell;
        torque = cross(displacement_vector1, forceV);
        HipMoment_1 = [-1*(Data(tt).LHipMoment(:, 1)), Data(tt).LHipMoment(:, 2), Data(tt).LHipMoment(:, 3)];
        torque_truncated = torque(1:size(HipMoment_1, 1), :);
        Total_torque= torque_truncated + HipMoment_1;
        hippos=Total_torque(:, 1);
        %plot(hipos');
        hippos1 = 1*(Data(tt).LGRF(:, 1));
        hippos = filtfilt(b, a, hippos);
    end
    
    % Use the pre-defined onset from the onsets variable
    onset_value = onsets(tt);
    
    %{
    % Plot the trial
    figure;
    plot(hippos1);
    xlabel('GRF');
    ylabel('Value');
    title(['Trial ', num2str(tt)]);
    
    % Mark the onset on the plot
    hold on;
    plot(onset_value, hippos1(onset_value), 'ro'); % Mark the onset with a red circle
    hold off;
    
    % Pause to visually confirm the onset (optional)
    pause(1); % Adjust the pause duration as needed
    close; % Close the figure after displaying
    %}
    Data(tt).HipAnkle1 = hippos;
end
% Subtract the mean of the first 10 data points from each column
hipposs = zeros(size(Data, 2), 571); % Initialize hipposs matrix
for tt = 1:size(Data, 2)
    means = mean(Data(tt).HipAnkle1(1:100));
    Data(tt).HipAnkle1 = Data(tt).HipAnkle1 - repmat(means, size(Data(tt).HipAnkle1, 1), 1);
    hipposs(tt,:) =   Data(tt).HipAnkle1(onsets(tt)-200:onsets(tt)+370); % based on onset take -75 and +500 points
    %hippos_stance_fro(tt,:)=hippos_stance_fro(onsets(tt)-170:onsets(tt)+500);
    %hippos_stance_sag(tt,:)=hippos_stance_sag(onsets(tt)-170:onsets(tt)+500);
end
%figure;
%plot(hipposs')

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

%{
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
%}
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

    %figure % Plotting onset
    %plot(hipacc)
    %hold on
    %plot([onset,onset],[min(hipacc),max(hipacc)],'b--')
end


% Initialize array to store onset values
load('a_ecu.mat', 'onsets1');

% Compute means of the first 10 data points for each trial
for tt = 1:size(Data, 2)
    if leg_choose == 1
        loadcell1 = 1*(Data(tt).RLoadCell);
        loadcell = (4.9984*(loadcell1) + 0.0062595)*4.4482216153;
        
        RB1_1 = [Data(tt).RB1(:, 1), Data(tt).RB1(:, 2), Data(tt).RB1(:, 3)];
        RB2_1 = [Data(tt).RB2(:, 1), Data(tt).RB2(:, 2), Data(tt).RB2(:, 3)];
        RFT_1= [Data(tt).RFT(:, 1), Data(tt).RFT(:, 2), Data(tt).RFT(:, 3)];
        
        % Calculate the displacement vector from RB1 to RB2
        displacement_vector = RB1_1 - RB2_1;
        displacement_vector1 = RFT_1;%- RB1_1;
        forceV= displacement_vector .* loadcell;
        torque = cross(displacement_vector1, forceV);
        HipMoment_1 = [Data(tt).RHipMoment(:, 1), Data(tt).RHipMoment(:, 2), Data(tt).RHipMoment(:, 3)];
        torque_truncated = torque(1:size(HipMoment_1, 1), :);
        Total_torque= torque_truncated + HipMoment_1;
        hipos=Total_torque(:, 1);
        %plot(hipos');
        hipos1 = Data(tt).RGRF(:, 1);

        %hipos = (4.9984*(hipos) + 0.0062595)*4.4482216153;
        %hippos_stance_fro=Data(tt).RHipAngle(:, 2);
        %hippos_stance_sag=filtfilt(b, a, hippos_stance_sag);
        %hippos_stance_fro=filtfilt(b, a, hippos_stance_fro);
        hipos = filtfilt(b, a, hipos);
    elseif leg_choose == 2
        loadcell1 = 1*(Data(tt).LLoadCell);
        loadcell = (5.0008*(loadcell1) + 0.0028211)*4.44822161533;
        
        RB1_1 = [(Data(tt).LB1(:, 1)), (Data(tt).LB1(:, 2)), (Data(tt).LB1(:, 3))];
        RB2_1 = [(Data(tt).LB2(:, 1)), (Data(tt).LB2(:, 2)), (Data(tt).LB2(:, 3))];
        RFT_1= [Data(tt).LFT(:, 1), Data(tt).LFT(:, 2), Data(tt).LFT(:, 3)];
        
        % Calculate the displacement vector from RB1 to RB2
        displacement_vector = RB1_1 - RB2_1;
        displacement_vector1 = RFT_1;%- RB1_1;
        forceV= displacement_vector .* loadcell;
        torque = cross(displacement_vector1, forceV);
        HipMoment_1 = [-1*(Data(tt).LHipMoment(:, 1)), Data(tt).LHipMoment(:, 2), Data(tt).LHipMoment(:, 3)];
        torque_truncated = torque(1:size(HipMoment_1, 1), :);
        Total_torque= torque_truncated + HipMoment_1;
        hipos=Total_torque(:, 1);
        %plot(hipos');
        hipos1 = 1*(Data(tt).LGRF(:, 1));
        hipos = filtfilt(b, a, hipos);
    end
    % Use the pre-defined onset from the onsets variable
    onset_value1 = onsets1(tt);
    %{
    % Plot the trial
    figure;
    plot(hipos1);
    xlabel('GRF');
    ylabel('Value');
    title(['Trial ', num2str(tt)]);
    
    % Mark the onset on the plot
    hold on;
    plot(onset_value1, hipos1(onset_value1), 'ro'); % Mark the onset with a red circle
    hold off;
    
    % Pause to visually confirm the onset (optional)
    pause(1); % Adjust the pause duration as needed
    close; % Close the figure after displaying
    %}
    Data(tt).HipAnkle11 = hipos;
end
% Subtract the mean of the first 10 data points from each column
hiposs = zeros(size(Data, 2), 571); % Initialize hipposs matrix
for tt = 1:size(Data, 2)
    means = mean(Data(tt).HipAnkle11(1:100));
    Data(tt).HipAnkle11 = Data(tt).HipAnkle11 - repmat(means, size(Data(tt).HipAnkle11, 1), 1);
    hiposs(tt,:) =   Data(tt).HipAnkle11(onsets1(tt)-200:onsets1(tt)+370); % based on onset take -75 and +500 points
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
%figure;
subplot(2,2,2)
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
%ylabel('Adduction (-) / Abduction (+) Torques (N.m)');
%title('1');
xticks(100:100:700);
xticklabels({'1','2','3','4','5','6','7'});
yticks(-80:20:40);
ylim([-80, 40]);
fontsize(16,"points");
hold off

% Save the variables into a MAT file
%save(fullfile(f0]older_path, 'max_noexo_right_sagittal.mat'), 'hipposs', 'm_hipposs', 'std_hipposs', 'onsets');


%folder_path = '/Users/mounikapasavula/Documents/Gait study /S01/StanceAngles';

% Save the variables into a MAT file in the specified folder
%save(fullfile(folder_path, 'max_normal_noexo_right_Offset.mat'), 'm_hipposs', 'std_hipposs', 'onsets', 'm_hiposs', 'std_hiposs', 'onsets1');

%% Repeating for third time

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

%{
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
%}
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

    %figure % Plotting onset
    %plot(hipacc)
    %hold on
    %plot([onset,onset],[min(hipacc),max(hipacc)],'b--')
end


% Initialize array to store onset values

load('a_noecu.mat', 'onsets');

% Compute means of the first 10 data points for each trial
for tt = 1:size(Data, 2)
    if leg_choose == 1
        hippos = 1 * (Data(tt).RHipMoment(:, 1));
        hippos1 = Data(tt).RGRF(:, 1);
        hippos = filtfilt(b, a, hippos);
    elseif leg_choose == 2
        hippos = -1* (Data(tt).RHipMoment(:, 1));
        %hippos=(5.0008*(hippos_) + 0.0028211)*4.44822161533;
        hippos1 = Data(tt).LGRF(:, 1);
        hippos = filtfilt(b, a, hippos);
    end
    
    % Use the pre-defined onset from the onsets variable
    onset_value = onsets(tt);

    %{
    % Plot the trial
    figure;
    plot(hippos1);
    xlabel('GRF');
    ylabel('Value');
    title(['Trial ', num2str(tt)]);
    
    % Mark the onset on the plot
    hold on;
    plot(onset_value, hippos1(onset_value), 'ro'); % Mark the onset with a red circle
    hold off;
    
    % Pause to visually confirm the onset (optional)
    pause(1); % Adjust the pause duration as needed
    close; % Close the figure after displaying
    %}
    Data(tt).HipAnkle1 = hippos;
end
% Subtract the mean of the first 10 data points from each column
hipposs = zeros(size(Data, 2), 571); % Initialize hipposs matrix
for tt = 1:size(Data, 2)
    means = mean(Data(tt).HipAnkle1(1:100));
    Data(tt).HipAnkle1 = Data(tt).HipAnkle1 - repmat(means, size(Data(tt).HipAnkle1, 1), 1);
    hipposs(tt,:) =   Data(tt).HipAnkle1(onsets(tt)-200:onsets(tt)+370); % based on onset take -75 and +500 points
    %hippos_stance_fro(tt,:)=hippos_stance_fro(onsets(tt)-170:onsets(tt)+500);
    %hippos_stance_sag(tt,:)=hippos_stance_sag(onsets(tt)-170:onsets(tt)+500);
end
%figure;
%plot(hipposs')

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

%{
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
%}
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

    %figure % Plotting onset
    %plot(hipacc)
    %hold on
    %plot([onset,onset],[min(hipacc),max(hipacc)],'b--')
end


% Initialize array to store onset values
load('a_noecu.mat', 'onsets1');

% Compute means of the first 10 data points for each trial
for tt = 1:size(Data, 2)
    if leg_choose == 1
        hipos = 1*(Data(tt).RHipMoment(:, 2));
        hipos1 = Data(tt).RGRF(:, 1);
        %hipos = (4.9984*(hipos) + 0.0062595)*4.4482216153;
        %hippos_stance_fro=Data(tt).RHipAngle(:, 2);
        %hippos_stance_sag=filtfilt(b, a, hippos_stance_sag);
        %hippos_stance_fro=filtfilt(b, a, hippos_stance_fro);
        hipos = filtfilt(b, a, hipos);
    elseif leg_choose == 2
        hipos = -1*(Data(tt).RHipMoment(:, 1));
        %hipos= (5.0008*(hipos_) + 0.0028211)*4.44822161533;
        hipos1 = Data(tt).LGRF(:, 1);
        %hipos = (5.0008*(hipos) + 0.0028211)*4.4482216153;
        %hippos_stance_sag=Data(tt).LHipAnkle(:, 1);
        %hippos_stance_fro=Data(tt).LHipAnkle(:, 2);
        %%hippos_stance_fro=filtfilt(b, a, hippos_stance_fro);
        hipos = filtfilt(b, a, hipos);
    end
    % Use the pre-defined onset from the onsets variable
    onset_value1 = onsets1(tt);

    %{
    % Plot the trial
    figure;
    plot(hipos1);
    xlabel('GRF');
    ylabel('Value');
    title(['Trial ', num2str(tt)]);
    
    % Mark the onset on the plot
    hold on;
    plot(onset_value1, hipos1(onset_value1), 'ro'); % Mark the onset with a red circle
    hold off;
    
    % Pause to visually confirm the onset (optional)
    pause(1); % Adjust the pause duration as needed
    close; % Close the figure after displaying
    %}
    Data(tt).HipAnkle11 = hipos;
end
% Subtract the mean of the first 10 data points from each column
hiposs = zeros(size(Data, 2), 571); % Initialize hipposs matrix
for tt = 1:size(Data, 2)
    means = mean(Data(tt).HipAnkle11(1:100));
    Data(tt).HipAnkle11 = Data(tt).HipAnkle11 - repmat(means, size(Data(tt).HipAnkle11, 1), 1);
    hiposs(tt,:) =   Data(tt).HipAnkle11(onsets1(tt)-200:onsets1(tt)+370); % based on onset take -75 and +500 points
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
%figure;
subplot(2,2,3)
plot(m_hipposs, 'color',[83/255, 27/255, 147/255] ,'LineWidth', 3.5);
%plot(m_hipposs, 'color',[148/255, 17/255, 0] ,'LineWidth', 3.5);
hold on;
x = 1:numel(m_hipposs);
fill([x fliplr(x)], [m_hipposs - std_hipposs fliplr(m_hipposs + std_hipposs)], [83/255, 27/255, 147/255], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
%fill([x fliplr(x)], [m_hipposs - std_hipposs fliplr(m_hipposs + std_hipposs)], [148/255, 17/255, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);


% Plot the second plot with its standard deviation shaded in lighter orange

plot(m_hiposs, 'Color', [0, 84/255, 147/255], 'LineWidth', 3.5); % Orange color
%plot(m_hiposs, 'Color', [255/255, 147/255, 0], 'LineWidth', 3.5);
hold on;
x = 1:numel(m_hiposs);
fill([x fliplr(x)], [m_hiposs - std_hiposs fliplr(m_hiposs + std_hiposs)], [0, 84/255, 147/255], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
%fill([x fliplr(x)], [m_hiposs - std_hiposs fliplr(m_hiposs + std_hiposs)], [255/255, 147/255, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
xlabel('Time (s)');
%ylabel('Adduction (-) / Abduction (+) Angles ( °)');
%ylabel('Assistance (-) / Resistance (+) Forces (N)');
%ylabel('Adduction (-) / Abduction (+) Torques (N.m)');
%title('1');
xticks(100:100:700);
xticklabels({'1','2','3','4','5','6','7'});
yticks(-60:20:80);
ylim([-60, 90]);
hold off
% Save the variables into a MAT file
%save(fullfile(f0]older_path, 'max_noexo_right_sagittal.mat'), 'hipposs', 'm_hipposs', 'std_hipposs', 'onsets');


%folder_path = '/Users/mounikapasavula/Documents/Gait study /S01/StanceAngles';

% Save the variables into a MAT file in the specified folder
%save(fullfile(folder_path, 'max_normal_noexo_right_Offset.mat'), 'm_hipposs', 'std_hipposs', 'onsets', 'm_hiposs', 'std_hiposs', 'onsets1');

%% Repeating for 4th time
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
%{
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
%}
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

    %figure % Plotting onset
    %plot(hipacc)
    %hold on
    %plot([onset,onset],[min(hipacc),max(hipacc)],'b--')
end


% Initialize array to store onset values

load('a_ecu.mat', 'onsets');

% Compute means of the first 10 data points for each trial
for tt = 1:size(Data, 2)
    if leg_choose == 1
        hippos = 1 * (Data(tt).RHipMoment(:, 1));
        hippos1 = Data(tt).RGRF(:, 1);
        hippos = filtfilt(b, a, hippos);
    elseif leg_choose == 2
        hippos = -1* (Data(tt).RHipMoment(:, 1));
        %hippos=(5.0008*(hippos_) + 0.0028211)*4.44822161533;
        hippos1 = Data(tt).LGRF(:, 1);
        hippos = filtfilt(b, a, hippos);
    end
    
    % Use the pre-defined onset from the onsets variable
    onset_value = onsets(tt);
%{
    % Plot the trial
    figure;
    plot(hippos1);
    xlabel('GRF');
    ylabel('Value');
    title(['Trial ', num2str(tt)]);
    
    % Mark the onset on the plot
    hold on;
    plot(onset_value, hippos1(onset_value), 'ro'); % Mark the onset with a red circle
    hold off;
    
    % Pause to visually confirm the onset (optional)
    pause(1); % Adjust the pause duration as needed
    close; % Close the figure after displaying
%}
    Data(tt).HipAnkle1 = hippos;
end
% Subtract the mean of the first 10 data points from each column
hipposs = zeros(size(Data, 2), 571); % Initialize hipposs matrix
for tt = 1:size(Data, 2)
    means = mean(Data(tt).HipAnkle1(1:100));
    Data(tt).HipAnkle1 = Data(tt).HipAnkle1 - repmat(means, size(Data(tt).HipAnkle1, 1), 1);
    hipposs(tt,:) =   Data(tt).HipAnkle1(onsets(tt)-200:onsets(tt)+370); % based on onset take -75 and +500 points
    %hippos_stance_fro(tt,:)=hippos_stance_fro(onsets(tt)-170:onsets(tt)+500);
    %hippos_stance_sag(tt,:)=hippos_stance_sag(onsets(tt)-170:onsets(tt)+500);
end
%figure;
%plot(hipposs')

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

%{
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
%}
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

    %figure % Plotting onset
    %plot(hipacc)
    %hold on
    %plot([onset,onset],[min(hipacc),max(hipacc)],'b--')
end


% Initialize array to store onset values
load('a_ecu.mat', 'onsets1');

% Compute means of the first 10 data points for each trial
for tt = 1:size(Data, 2)
    if leg_choose == 1
        hipos = 1*(Data(tt).RHipMoment(:, 2));
        hipos1 = Data(tt).RGRF(:, 1);
        %hipos = (4.9984*(hipos) + 0.0062595)*4.4482216153;
        %hippos_stance_fro=Data(tt).RHipAngle(:, 2);
        %hippos_stance_sag=filtfilt(b, a, hippos_stance_sag);
        %hippos_stance_fro=filtfilt(b, a, hippos_stance_fro);
        hipos = filtfilt(b, a, hipos);
    elseif leg_choose == 2
        hipos = -1*(Data(tt).RHipMoment(:, 1));
        %hipos= (5.0008*(hipos_) + 0.0028211)*4.44822161533;
        hipos1 = Data(tt).LGRF(:, 1);
        %hipos = (5.0008*(hipos) + 0.0028211)*4.4482216153;
        %hippos_stance_sag=Data(tt).LHipAnkle(:, 1);
        %hippos_stance_fro=Data(tt).LHipAnkle(:, 2);
        %%hippos_stance_fro=filtfilt(b, a, hippos_stance_fro);
        hipos = filtfilt(b, a, hipos);
    end
    % Use the pre-defined onset from the onsets variable
    onset_value1 = onsets1(tt);

    %{
    % Plot the trial
    figure;
    plot(hipos1);
    xlabel('GRF');
    ylabel('Value');
    title(['Trial ', num2str(tt)]);
    
    % Mark the onset on the plot
    hold on;
    plot(onset_value1, hipos1(onset_value1), 'ro'); % Mark the onset with a red circle
    hold off;
    
    % Pause to visually confirm the onset (optional)
    pause(1); % Adjust the pause duration as needed
    close; % Close the figure after displaying
    %}
    Data(tt).HipAnkle11 = hipos;
end
% Subtract the mean of the first 10 data points from each column
hiposs = zeros(size(Data, 2), 571); % Initialize hipposs matrix
for tt = 1:size(Data, 2)
    means = mean(Data(tt).HipAnkle11(1:100));
    Data(tt).HipAnkle11 = Data(tt).HipAnkle11 - repmat(means, size(Data(tt).HipAnkle11, 1), 1);
    hiposs(tt,:) =   Data(tt).HipAnkle11(onsets1(tt)-200:onsets1(tt)+370); % based on onset take -75 and +500 points
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

subplot(2,2,4)
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
%ylabel('Adduction (-) / Abduction (+) Torques (N.m)');
%title('1');
xticks(100:100:700);
xticklabels({'1','2','3','4','5','6','7'});
yticks(-60:20:80);
ylim([-60, 90]);
fontsize(16,"points");
hold off
% Save the variables into a MAT file
%save(fullfile(f0]older_path, 'max_noexo_right_sagittal.mat'), 'hipposs', 'm_hipposs', 'std_hipposs', 'onsets');


%folder_path = '/Users/mounikapasavula/Documents/Gait study /S01/StanceAngles';

% Save the variables into a MAT file in the specified folder
%save(fullfile(folder_path, 'max_normal_noexo_right_Offset.mat'), 'm_hipposs', 'std_hipposs', 'onsets', 'm_hiposs', 'std_hiposs', 'onsets1');
