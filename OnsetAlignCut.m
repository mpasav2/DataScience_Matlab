clear; clc;
%% Onset detection
load S02_cut_newonset.mat % combined and manually fixed version
% offset
trialnums = [6 7 6 7 7 6 7 7];
trial = 1;
baserange = 1:150;
detectrange = 1:500;

onset = [];
for cc=1:8
    for tt=1:trialnums(cc)
        if cc == 1 || cc == 2 || cc == 5 || cc == 6
            hipflex = S02(trial).LHipAngle(1,:);   % flexion of the stepping leg
        elseif cc == 3 || cc == 4 || cc == 7 || cc == 8
            hipflex = S02(trial).RHipAngle(1,:);
        end

        % offset by max extension angle in stance leg in each trial
        onset(trial,1) = find(hipflex==min(hipflex));
%         % offset by mean of baseline in each trial
%         grfoff = [S02(trial).RGRF;S02(trial).LGRF] - repmat(mean([S02(trial).RGRF(:,baserange);S02(trial).LGRF(:,baserange)],2),1,size(S02(trial).RGRF,2));
%         % composite metric which detects the switching of the weight shift
%         % during automatic postural adjustment
%         grfcomp = abs(grfoff(3,detectrange)).*abs(grfoff(6,detectrange)); 
%         [pks,locs] = findpeaks(grfcomp,'MinPeakHeight',0.005,'MinPeakWidth',10);    % finding peaks

        % either the point at which the first peak from above or at 202, to
        % capture at least 2s of standing, minus the very first two samples
%         onset(trial,1) = max(locs(1),202);  % turns out only two trials, each 1 or 2 samples away => So, we're OK

        % Visually inspect
        subplot(8,7,7*(cc-1)+tt)
%         plot(detectrange,grfcomp)
%         hold on
%         plot(locs, pks, 'x')
%         yyaxis right
        plot(detectrange,hipflex(detectrange))
        hold on
        plot([onset(trial) onset(trial)],ylim,'k')
        box off;

        if tt==1
            ylabel([S02(trial).Exo '-' S02(trial).Leg '-' S02(trial).Step])
        end
        trial = trial+1;
    end
end

%% Filter, resample, convert, and cut
clc;

Fs_emg = 290/0.135; % Fs(1); % emg sampling frequency, weird number due to sampling settings; can load from one of the recorded files 
Fs_mot = 100; % motion capture sampling frequency

emgtime = ((1:24000)-1)/Fs_emg;   % based on longest trial
mottime = ((1:1501)-1)/Fs_mot;  % based on longest trial

% low-pass filter
F_low = 10; % cutoff
[bemg,aemg] = butter(2, F_low/(Fs_emg/2), 'low');
[bmot,amot] = butter(2, F_low/(Fs_mot/2), 'low');

MotVars = {'RHipAngle','RKneeAngle','RAnkleAngle', ... % right to left, angles, 
             'LHipAngle','LKneeAngle','LAnkleAngle', ...
             'ThoraxAngle', ...
             'RHipMoment','RKneeMoment','RAnkleMoment', ... % Moments
             'LHipMoment','LKneeMoment','LAnkleMoment', ...
             'RGRF', 'LGRF', 'RLoadCell', 'LLoadCell'};     % External forces

for tt=1:size(S02,2)
    onset_trial = onset(tt);
    ind_trial = (onset_trial-249):(onset_trial+350);
    time_trial = mottime(ind_trial);

    % Conditions
    S02cut(tt).Exo = S02(tt).Exo;
    S02cut(tt).Leg = S02(tt).Leg;
    S02cut(tt).Step = S02(tt).Step;

    % EMGs
    remg = S02(tt).REMG;
    remgfilt = filtfilt(bemg,aemg,abs(remg)')';
    for mm=1:8
        S02cut(tt).REMG(mm,:) = interp1(emgtime(1:length(remg(mm,:))), remgfilt(mm,:), time_trial,'pchip');
    end
    lemg = S02(tt).LEMG;
    lemgfilt = filtfilt(bemg,aemg,abs(lemg)')';
    for mm=1:8
        S02cut(tt).LEMG(mm,:) = interp1(emgtime(1:length(lemg(mm,:))), lemgfilt(mm,:), time_trial,'pchip');
    end

    % kinematic and kinetic variables
    for vv=1:15
        eval(['kvar = S02(tt).' MotVars{vv} ';'])
        if vv==7 % fixing trunk angle
            kvar(1,kvar(1,:)<0) = 360 + kvar(1,kvar(1,:)<0);
        elseif vv>7 && vv<14 % filtering just the moments
            kvar = filtfilt(bmot,amot,kvar')';
        end
        eval(['S02cut(tt).' MotVars{vv} ' = kvar(:,ind_trial);'])
    end

    if tt>26
        % % convert loadcell data to forces in N, and filter        % 
        % % calibration table  for right (dotted) loadcell: CertNo 2310230055
        % %  Load (lb) Output (Vdc)
        % calR = [0, 0; ...
        %         10, 1.998; ...
        %         20, 3.995; ...
        %         30, 5.999; ...
        %         40, 8.003; ...
        %         50, 10.004; ...
        %         0, 0.002];        % 
        % 
        % % for left loadcell: CertNo 2310230056
        % calㅣ = [0, 0; ...
        %         10, 1.998; ...
        %         20, 3.994; ...
        %         30, 6.000; ...
        %         40, 7.998; ...
        %         50, 9.999; ...
        %         0, 0.001];
        % 
        % % Could to linear regression. But coefficient can be simply considered to be 0.2
        % 
        scale = 0.2*4.4482216153; % lbf to N conversion factor
        
        % loadcell data
        rload = filtfilt(bmot,amot,S02(tt).RLoadCell);
        llaod = filtfilt(bmot,amot,S02(tt).LLoadCell);
        S02cut(tt).RLoadCell = scale*rload(:,ind_trial);
        S02cut(tt).LLoadCell = scale*llaod(:,ind_trial);
    end
end

% saved to S02_cut.mat





