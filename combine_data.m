% clear; clc;
% list []condition labelss to reorganize (combine data)
step = {'max step', 'normal step'};
leg = {'left', 'right'};
trialnum = [5, 4, 5, 4];
cumnum = cumsum(trialnum);

%datafolder = 'S04'; % folder contatining trial data
%cd(datafolder)

trial = 1;
for ss=1:length(step)
    for ll=1:length(leg)
        for tt=1:trialnum(length(step)*(ss-1)+ll)
            filename = [step{ss} ' ' leg{ll} '_' int2str(tt) '.mat'];
            load(filename)
%             if ss==1 || ll==1 || tt==1
% %                 S = whos;
% %                 varnames = [];
% %                 for vv=1:22 %size(S,1)
% %                     varnames{vv} = S(vv).name;
% %                 end
%                 % variable selection in order of preference
                vars = [9, 11, 6, 16, 18, 14, 22, ... % joint kinematics: leg (left to right, proximal to distal), then thorax
                        10, 13, 7, 17, 19, 15, ... % joing moments
                        2, 3, ... COM COP
                        12, 20]; %
%                 % clear S
%             end
            if ss==1 && ll==1
                trial = tt;
            else
                trial = cumnum(length(step)*(ss-1)+ll-1) + tt;
            end
            Data(trial).Cond = [step{ss}(1:(length(step{ss})-5)) '_' leg{ll} int2str(tt)];
            for vv = 1:length(vars)
                eval(['Data(trial).' varnames{vars(vv)} '=  ' varnames{vars(vv)} '{1,1};']);
            end
        end
    end
end
cd .. % move back to parent folder

Fs_mot = 100; % motion capture sampling frequency
Fs_load = 1000; % loadcell sampling frequency


%% converting signs, based on axis (in mocap) vs.  
for tt=1:18
    Data(tt).LHipAngle(:,2) = -Data(tt).LHipAngle(:,2);     % hip ab/adduction in the left leg
    Data(tt).LHipMoment(:,2) = -Data(tt).LHipMoment(:,2);
    
    Data(tt).LKneeAngle(:,2) = -Data(tt).LKneeAngle(:,2); % knee frontal plane angle in the left leg
    Data(tt).LKneeMoment(:,2) = -Data(tt).LKneeMoment(:,2);
    
    Data(tt).LAnkleAngle(:,2) = -Data(tt).LAnkleAngle(:,2); % Ankle frontal plane angle in the left leg
    Data(tt).LAnkleMoment(:,2) = -Data(tt).LAnkleMoment(:,2);
end

%% plot

% frontal plan hip flexion angles
% to check which variable to use to detect the onset of movement
close all; clc;
for tt=1:5
    data = Data(tt).LHipAngle(:,1);
    time = (1:length(data))/Fs_mot;

    subplot(3,2,tt)
    yyaxis left
    plot(time, data)
    yyaxis right
    plot(time, gradient(gradient(data,time),time));
end

figure;
for tt=6:9
    data = Data(tt).RHipAngle(:,1);
    time = (1:length(data))/Fs_mot;

    subplot(3,2,tt-5)
    yyaxis left
    plot(time, data)
    yyaxis right
    plot(time, gradient(gradient(data,time),time));
end

figure;
for tt=10:14
    data = Data(tt).LHipAngle(:,1);
    time = (1:length(data))/Fs_mot;

    subplot(3,2,tt-9)
    yyaxis left
    plot(time, data)
    yyaxis right
    plot(time, gradient(gradient(data,time),time));
end

figure;
for tt=15:18
    data = Data(tt).RHipAngle(:,1);
    time = (1:length(data))/Fs_mot;

    subplot(3,2,tt-14)
    yyaxis left
    plot(time, data)
    yyaxis right
    plot(time, gradient(gradient(data,time),time));
end

% We conclude that using first inflection point of hip flexion acc could
% work
%% Onset detection

% Low-pass filtering
F_low = 4; % cutoff lower bound
[b,a] = butter(2, F_low/(Fs_mot/2), 'low');

for tt=1:18
    if strcmp(Data(tt).Cond((end-2)),'f')
        hippos = Data(tt).LHipAngle(:,1);
    elseif strcmp(Data(tt).Cond((end-2)),'h')
        hippos = Data(tt).RHipAngle(:,1);
    end
    time = (1:length(hippos))/Fs_mot;

    hipvel = gradient(hippos,time);
    hipvel = filtfilt(b,a,hipvel);

    hipacc = gradient(hipvel,time);
    hipacc = filtfilt(b,a,hipacc);
    
    Data(tt).hipvel = hipvel;
    Data(tt).hipacc = hipacc;
end

close all; clc;
for tt=1:5
    data = Data(tt).LHipAngle(:,1);
    time = (1:length(data))/Fs_mot;

    subplot(3,2,tt)
    yyaxis left
    plot(time, data)
    yyaxis right
    plot(time, Data(tt).hipacc);
end

figure;
for tt=6:9
    data = Data(tt).RHipAngle(:,1);
    time = (1:length(data))/Fs_mot;

    subplot(3,2,tt-5)
    yyaxis left
    plot(time, data)
    yyaxis right
    plot(time, Data(tt).hipacc);
end

figure;
for tt=10:14
    data = Data(tt).LHipAngle(:,1);
    time = (1:length(data))/Fs_mot;

    subplot(3,2,tt-9)
    yyaxis left
    plot(time, data)
    yyaxis right
    plot(time, Data(tt).hipacc);
end

figure;
for tt=15:18
    data = Data(tt).RHipAngle(:,1);
    time = (1:length(data))/Fs_mot;

    subplot(3,2,tt-14)
    yyaxis left
    plot(time, data)
    yyaxis right
    plot(time, Data(tt).hipacc);
end






