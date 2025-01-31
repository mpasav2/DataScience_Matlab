clear; clc; close all;
%load S02_cut
load S02_cut_newonset

musnames = {'Gmax','Gmed','Adl','RF','VL','BF','TA','MG'}; % muscle names
% conditions
exo = {'WO', 'W'};
leg = {'R', 'L'};
step = {'Norm', 'Long'};

trialnums = [6 7 6 7 7 6 7 7]; % number of trials/repetitions for each condition

% time = (-199:400)/Fs_mot; % onset as t=0, then -2s to +4s (6s in total)
time = (-249:350)/Fs_mot; 
%% Plot for each condition
%close all;
figure;

% Choose options to plot
var = 'REMG(1,:)'; % choose the variable; also specify muscle or DoF
% examples: var = 'LEMG(2,:)' for Gmed muscle of Left leg
% examples: var = 'RHipAngle(1,:)' for flexion DoF of Right leg
opt = 'avg';    % choose whether to plot individual trials or avg+/-std
baselinecorrect = 0; % choose whether to remove the baseline (first 1.5s) or not
datarange = 1:400; % choose, in index, what range of data to plot, 1:600 is full

output = 0; % choose to output (if 1, not if 0) the variable for each condition
print2pdf = 0; % choose to print to a pdf file, filename as the variable name; will also close all figures not to be confused



% Initialize
datamin = 0;    % 0 doesn't work if all data points are >0
datamax = 0;

if strcmp(var(2:5),'Load')
    exorange = 2:2; % Loadcell data only in W Exo conditions
    trial = 27;
else
    exorange = 1:length(exo);
    trial=1;
end
% Output variables
output_data = struct();
for ee=exorange
    for ll = 1:length(leg)
        for ss = 1:length(step)
            cond = length(leg)*length(step)*(ee-1) + length(step)*(ll-1) + ss;
            data = [];
            subplot(4,2,cond)
            switch opt
                case 'ind'
                    for tt=1:trialnums(cond)
                        eval(['data = S02cut(trial).' var ';'])
                        if baselinecorrect
                            data = data - mean(data(1:150));    % up to 1.5s
                        end
                        data = data(:,datarange);
                        if ee*ll*ss*tt == 1
                            if min(data) > 0
                                datamin = min(data);
                            end
                        end

                        if output
                            eval([var(1:7) '_' exo{ee} '_' leg{ll} '_' step{ss} int2str(tt) ' = data;'])
                        end

                        plot(time(datarange),data,'Color',[0.5 0.5 0.5])
                        hold on
                        trial = trial+1;
                        datamin = min(datamin, min(data));
                        datamax = max(datamax, max(data));
                    end
                case 'avg'
                    for tt=1:trialnums(cond)
                        eval(['data = [data; S02cut(trial).' var '];'])
                        trial = trial+1;
                    end
                    if baselinecorrect
                        data = data - repmat(mean(data(:,1:150),2),1,size(data,2));
                    end
                    data = data(:,datarange);
                    if ee*ll*ss == 1
                        if min(min(data)) > 0
                            datamin = min(min(data));
                        end
                    end

                    if output
                        eval([var(1:7) '_' exo{ee} '_' leg{ll} '_' step{ss} ' = data;'])
                    end

                    fill([time(datarange), fliplr(time(datarange))], ...
                        [mean(data,1)+std(data,0,1), fliplr(mean(data,1)-std(data,0,1))], ...
                        [0.75 0.75 0.75],'EdgeColor','none')
                    hold on
                    plot(time(datarange),mean(data,1),'k','LineWidth',1.25)
                    datamin = min(datamin, min(mean(data,1)-std(data,1,1)));
                    datamax = max(datamax, max(mean(data,1)+std(data,1,1)));
                    mean_data = mean(data(:, datarange), 1);
                    std_data = std(data(:, datarange), 0, 1);
                    output_data(cond).mean = mean_data;
                    output_data(cond).std = std_data;
                    output_data(cond).onset = onset(cond);
            end
            box off; 
            if ee==1 && ll==1
                title(step{ss})
            end
            if ss==1
                if ee==1
                    ylabel(['w/o Exo --- Step w/' S02cut(trial).Leg])
                elseif ee==2
                    ylabel(['w Exo --- Step w/' S02cut(trial).Leg])
                end
            end
        end
    end
end
% put all in the same scale
for cc=1:8
    subplot(4,2,cc)
    %ylim([170, 190])
    if datamin < 0
        ylim([1.05*datamin, 1.05*datamax])
    else
        ylim([0.95*datamin, 1.05*datamax])
    end
    plot([0 0],[datamin datamax],':r')
end

if print2pdf
    print(var,'-dpdf','-fillpage')
    close all
end


%onset_data = load('S02_cut_newonset.mat');
%S02_cut_newonset = onset.S02_cut_newonset;
% Define folder path
%folder_path = '/Users/mounikapasavula/Desktop/Hongchul Code/New codes';

% Save the variables into a MAT file in the specified folder
%save(fullfile(folder_path, 'test.mat'), 'output_data','S02_cut_newonset');
