function S01 = EMG_import(params)

% Moving to the folder with the data
cd(params.EMG_datafolder);
% baseline trials
    trial = 1;
    trials = []; trials = params.trialnum{1};
    for tt=1:length(trials)
        filename = ['Run_number_35_Plot_and_Store_Rep_' int2str(trials(tt)) '.' int2str(trials(tt)-1) '.mat'];
        load(filename)
        S01(trial).Exo = 'base';
        S01(trial).Leg = 'base';
        S01(trial).Step = 'base';
        %indnonzero = find(Data(1,:)~=0 & Data(end,:)~=0); % first and last, just in case
        emg_r = []; emg_r = detrend(Data(params.muslist(:,1),:)')';
        emg_l = []; emg_l = detrend(Data(params.muslist(:,2),:)')';
        S01(trial).REMG = emg_r;
        S01(trial).LEMG = emg_l;
        S01(trial).time = Time(1,1:length(Data));
        trial = trial+1;
        clear Data Time Fs
    end
    
    for ee=1:length(params.exo)
        for ll=1:length(params.leg)
            for ss=1:length(params.step)
                cond = length(params.leg)*length(params.step)*(ee-1) + length(params.step)*(ll-1) + ss + 1;
                trials= params.trialnum{cond};
                for tt=1:length(trials)
                    filename = ['Run_number_35_Plot_and_Store_Rep_' int2str(trials(tt)) '.' int2str(trials(tt)-1) '.mat'];
                    load(filename)
                    S01(trial).Exo = params.exo{ee};
                    S01(trial).Leg = params.leg{ll};
                    S01(trial).Step = params.step{ss};
                    %indnonzero = find(Data(1,:)~=0 & Data(end,:)~=0); % first and last, just in case
                    emg_r = []; emg_r = detrend(Data(params.muslist(:,1),:)')';
                    emg_l = []; emg_l = detrend(Data(params.muslist(:,2),:)')';
                    S01(trial).REMG = emg_r;
                    S01(trial).LEMG = emg_l;
                    S01(trial).time = Time(1,1:length(Data));
                    trial = trial+1;
                    clear Data Time Fs
                end
            end
        end
    end

    cd .. % move back to parent folder
end
