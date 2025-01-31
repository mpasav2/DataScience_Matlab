% Define conditions
exo = {'BA', 'noexo', 'noexo', 'exo', 'exo'};
leg = {'N/A', 'left', 'left', 'left', 'left'};
step = {'N/A', 'normal', 'max', 'normal', 'max'};

% Define trial numbers for each condition
trialnum = {
    1:3, ... % Baseline trials
    4:10, ... % WO R Norm
    11:17, ... % W0 R Long
    18:24, ... % W R Norm
    25:31, ... % W R Long
};
trial_index = 1;
for i = 1:length(S02)

    if ismember(S02(i).Exo,exo{1})
        continue;
    else
            S02(i).RCh1{1,:} = filter_segment(S02(i).Right{1,:},S02(i).onset,200,500,S02(i).Time{1,:});
            S02(i).RCh2{1,:} = filter_segment(S02(i).Right{2,:},S02(i).onset,200,500,S02(i).Time{2,:});
            S02(i).RCh3{1,:} = filter_segment(S02(i).Right{3,:},S02(i).onset,200,500,S02(i).Time{3,:});
            S02(i).RCh4{1,:} = filter_segment(S02(i).Right{4,:},S02(i).onset,200,500,S02(i).Time{4,:});
            S02(i).RCh5{1,:} = filter_segment(S02(i).Right{5,:},S02(i).onset,200,500,S02(i).Time{5,:});
            S02(i).RCh6{1,:} = filter_segment(S02(i).Right{6,:},S02(i).onset,200,500,S02(i).Time{6,:});
            S02(i).RCh7{1,:} = filter_segment(S02(i).Right{7,:},S02(i).onset,200,500,S02(i).Time{7,:});
            S02(i).RCh8{1,:} = filter_segment(S02(i).Right{8,:},S02(i).onset,200,500,S02(i).Time{8,:});

            S02(i).LCh1{1,:} = filter_segment(S02(i).Left{1,:},S02(i).onset,200,500,S02(i).Time{1,:});
            S02(i).LCh2{1,:} = filter_segment(S02(i).Left{2,:},S02(i).onset,200,500,S02(i).Time{2,:});
            S02(i).LCh3{1,:} = filter_segment(S02(i).Left{3,:},S02(i).onset,200,500,S02(i).Time{3,:});
            S02(i).LCh4{1,:} = filter_segment(S02(i).Left{4,:},S02(i).onset,200,500,S02(i).Time{4,:});
            S02(i).LCh5{1,:} = filter_segment(S02(i).Left{5,:},S02(i).onset,200,500,S02(i).Time{5,:});
            S02(i).LCh6{1,:} = filter_segment(S02(i).Left{6,:},S02(i).onset,200,500,S02(i).Time{6,:});
            S02(i).LCh7{1,:} = filter_segment(S02(i).Left{7,:},S02(i).onset,200,500,S02(i).Time{7,:});
            S02(i).LCh8{1,:} = filter_segment(S02(i).Left{8,:},S02(i).onset,200,500,S02(i).Time{8,:});
    end
end


for i = 1:length(S02)
    if S02(i).onset ~= 0
        S02(i).GMax = filtering(S02(i).Right{1,1},S02(i).Hip_torque(:,1),S02(i).Time{1,1},1,0);
        S02(i).GMed = filtering(S02(i).Right{2,1},S02(i).Hip_torque(:,2),S02(i).Time{2,1},2,0);
        S02(i).AddL = filtering(S02(i).Right{3,1},S02(i).Hip_torque(:,2),S02(i).Time{3,1},1,0);
    end
end
%%
S03 = struct();

for i = 1:length(trialnum)
    if i ~= 1
        Chs_R1 = [];
        Chs_L1 = [];
        Chs_R2 = [];
        Chs_L2 = [];
        Chs_R3 = [];
        Chs_L3 = [];
        Chs_R4 = [];
        Chs_L4 = [];
        Chs_R5 = [];
        Chs_L5 = [];
        Chs_R6 = [];
        Chs_L6 = [];
        Chs_R7 = [];
        Chs_L7 = [];
        Chs_R8 = [];
        Chs_L8 = [];


        for j = 1:length(trialnum{i})
            indx_cond = trialnum{i};
            indx_trial = indx_cond(j);
            if size(S02(indx_trial).RCh1{1,1},1) == 1
                Chs_R1 = [Chs_R1; S02(indx_trial).RCh1{1,1}];
                Chs_L1 = [Chs_L1; S02(indx_trial).LCh1{1,1}];
    
                Chs_R2 = [Chs_R2; S02(indx_trial).RCh2{1,1}];
                Chs_L2 = [Chs_L2; S02(indx_trial).LCh2{1,1}];
    
                Chs_R3 = [Chs_R3; S02(indx_trial).RCh3{1,1}];
                Chs_L3 = [Chs_L3; S02(indx_trial).LCh3{1,1}];
    
                Chs_R4 = [Chs_R4; S02(indx_trial).RCh4{1,1}];
                Chs_L4 = [Chs_L4; S02(indx_trial).LCh4{1,1}];
    
                Chs_R5 = [Chs_R5; S02(indx_trial).RCh5{1,1}];
                Chs_L5 = [Chs_L5; S02(indx_trial).LCh5{1,1}];
    
                Chs_R6 = [Chs_R6; S02(indx_trial).RCh6{1,1}];
                Chs_L6 = [Chs_L6; S02(indx_trial).LCh6{1,1}];
    
                Chs_R7 = [Chs_R7; S02(indx_trial).RCh7{1,1}];
                Chs_L7 = [Chs_L7; S02(indx_trial).LCh7{1,1}];
    
                Chs_R8 = [Chs_R8; S02(indx_trial).RCh8{1,1}];
                Chs_L8 = [Chs_L8; S02(indx_trial).LCh8{1,1}];
            else
                Chs_R1 = [Chs_R1; S02(indx_trial).RCh1{1,1}'];
                Chs_L1 = [Chs_L1; S02(indx_trial).LCh1{1,1}'];
    
                Chs_R2 = [Chs_R2; S02(indx_trial).RCh2{1,1}'];
                Chs_L2 = [Chs_L2; S02(indx_trial).LCh2{1,1}'];
    
                Chs_R3 = [Chs_R3; S02(indx_trial).RCh3{1,1}'];
                Chs_L3 = [Chs_L3; S02(indx_trial).LCh3{1,1}'];
    
                Chs_R4 = [Chs_R4; S02(indx_trial).RCh4{1,1}'];
                Chs_L4 = [Chs_L4; S02(indx_trial).LCh4{1,1}'];
    
                Chs_R5 = [Chs_R5; S02(indx_trial).RCh5{1,1}'];
                Chs_L5 = [Chs_L5; S02(indx_trial).LCh5{1,1}'];
    
                Chs_R6 = [Chs_R6; S02(indx_trial).RCh6{1,1}'];
                Chs_L6 = [Chs_L6; S02(indx_trial).LCh6{1,1}'];
    
                Chs_R7 = [Chs_R7; S02(indx_trial).RCh7{1,1}'];
                Chs_L7 = [Chs_L7; S02(indx_trial).LCh7{1,1}'];
    
                Chs_R8 = [Chs_R8; S02(indx_trial).RCh8{1,1}'];
                Chs_L8 = [Chs_L8; S02(indx_trial).LCh8{1,1}'];
            end
        end
        Chs_R1avg = mean(Chs_R1,1);
        Chs_L1avg = mean(Chs_L1,1);

        Chs_R2avg = mean(Chs_R2,1);
        Chs_L2avg = mean(Chs_L2,1);

        Chs_R3avg = mean(Chs_R3,1);
        Chs_L3avg = mean(Chs_L3,1);

        Chs_R4avg = mean(Chs_R4,1);
        Chs_L4avg = mean(Chs_L4,1);

        Chs_R5avg = mean(Chs_R5,1);
        Chs_L5avg = mean(Chs_L5,1);

        Chs_R6avg = mean(Chs_R6,1);
        Chs_L6avg = mean(Chs_L6,1);

        Chs_R7avg = mean(Chs_R7,1);
        Chs_L7avg = mean(Chs_L7,1);

        Chs_R8avg = mean(Chs_R8,1);
        Chs_L8avg = mean(Chs_L8,1);

        Chs_R1std = std(Chs_R1,1);
        Chs_L1std = std(Chs_L1,1);

        Chs_R2std = std(Chs_R2,1);
        Chs_L2std = std(Chs_L2,1);

        Chs_R3std = std(Chs_R3,1);
        Chs_L3std = std(Chs_L3,1);

        Chs_R4std = std(Chs_R4,1);
        Chs_L4std = std(Chs_L4,1);

        Chs_R5std = std(Chs_R5,1);
        Chs_L5std = std(Chs_L5,1);

        Chs_R6std = std(Chs_R6,1);
        Chs_L6std = std(Chs_L6,1);

        Chs_R7std = std(Chs_R7,1);
        Chs_L7std = std(Chs_L7,1);

        Chs_R8std = std(Chs_R8,1);
        Chs_L8std = std(Chs_L8,1);

        S03(trial_index).Exo = exo{i};
        S03(trial_index).Leg = leg{i};
        S03(trial_index).Step = step{i};
        S03(trial_index).Chs_R1 = Chs_R1;
        S03(trial_index).Chs_L1 = Chs_L1;
        S03(trial_index).Chs_R2 = Chs_R2;
        S03(trial_index).Chs_L2 = Chs_L2;
        S03(trial_index).Chs_R3 = Chs_R3;
        S03(trial_index).Chs_L3 = Chs_L3;
        S03(trial_index).Chs_R4 = Chs_R4;
        S03(trial_index).Chs_L4 = Chs_L4;
        S03(trial_index).Chs_R5 = Chs_R5;
        S03(trial_index).Chs_L5 = Chs_L5;
        S03(trial_index).Chs_R6 = Chs_R6;
        S03(trial_index).Chs_L6 = Chs_L6;
        S03(trial_index).Chs_R7 = Chs_R7;
        S03(trial_index).Chs_L7 = Chs_L7;
        S03(trial_index).Chs_R8 = Chs_R8;
        S03(trial_index).Chs_L8 = Chs_L8;

        S03(trial_index).Ch_R1avg = Chs_R1avg;
        S03(trial_index).Ch_L1avg = Chs_L1avg;
        S03(trial_index).Ch_R2avg = Chs_R2avg;
        S03(trial_index).Ch_L2avg = Chs_L2avg;
        S03(trial_index).Ch_R3avg = Chs_R3avg;
        S03(trial_index).Ch_L3avg = Chs_L3avg;
        S03(trial_index).Ch_R4avg = Chs_R4avg;
        S03(trial_index).Ch_L4avg = Chs_L4avg;
        S03(trial_index).Ch_R5avg = Chs_R5avg;
        S03(trial_index).Ch_L5avg = Chs_L5avg;
        S03(trial_index).Ch_R6avg = Chs_R6avg;
        S03(trial_index).Ch_L6avg = Chs_L6avg;
        S03(trial_index).Ch_R7avg = Chs_R7avg;
        S03(trial_index).Ch_L7avg = Chs_L7avg;
        S03(trial_index).Ch_R8avg = Chs_R8avg;
        S03(trial_index).Ch_L8avg = Chs_L8avg;

        S03(trial_index).Ch_R1std = Chs_R1std;
        S03(trial_index).Ch_L1std = Chs_L1std;
        S03(trial_index).Ch_R2std = Chs_R2std;
        S03(trial_index).Ch_L2std = Chs_L2std;
        S03(trial_index).Ch_R3std = Chs_R3std;
        S03(trial_index).Ch_L3std = Chs_L3std;
        S03(trial_index).Ch_R4std = Chs_R4std;
        S03(trial_index).Ch_L4std = Chs_L4std;
        S03(trial_index).Ch_R5std = Chs_R5std;
        S03(trial_index).Ch_L5std = Chs_L5std;
        S03(trial_index).Ch_R6std = Chs_R6std;
        S03(trial_index).Ch_L6std = Chs_L6std;
        S03(trial_index).Ch_R7std = Chs_R7std;
        S03(trial_index).Ch_L7std = Chs_L7std;
        S03(trial_index).Ch_R8std = Chs_R8std;
        S03(trial_index).Ch_L8std = Chs_L8std;
        
        trial_index = trial_index + 1;
    end
end

%% Min max of plots

Chs_R_min = ones(8,1)*inf;
Chs_L_min = ones(8,1)*inf;
Chs_R_max = ones(8,1)*-inf;
Chs_L_max = ones(8,1)*-inf;

for n = 1:size(S03,2)
for i = 1:8
    if Chs_R_min(i) > min(min(S03(n).(['Chs_R' num2str(i)])))
        Chs_R_min(i) = min(min(S03(n).(['Chs_R' num2str(i)])));
    end
    if Chs_L_min(i) > min(min(S03(n).(['Chs_L' num2str(i)])))
        Chs_L_min(i) = min(min(S03(n).(['Chs_L' num2str(i)])));
    end
    if Chs_R_max(i) < max(max(S03(n).(['Chs_R' num2str(i)])))
        Chs_R_max(i) = max(max(S03(n).(['Chs_R' num2str(i)])));
    end
    if Chs_L_max(i) < max(max(S03(n).(['Chs_L' num2str(i)])))
        Chs_L_max(i) = max(max(S03(n).(['Chs_L' num2str(i)])));
    end
end
end
%% Plotting each trial

figure;
n = 1;
time = -2.7:1:5.0;

for i = 1:8
    subplot(4,2,i);
    time_Ch = linspace(-2.7, 5.0, length(S03(n).(['Chs_L' num2str(i)])));
    plot(time_Ch, S03(n).(['Chs_L' num2str(i)]), 'Color', [0.5 0.5 0.5]);
    hold on;
    plot(time_Ch, S03(n).(['Ch_L' num2str(i) 'avg']), 'k-', 'LineWidth', 2);
    ylim([Chs_L_min(i), Chs_L_max(i)]); % Set y-axis limits
    hold off;
end

%%
figure;
n = 1;
for i = 1:8
    subplot(4,2,i);
    time_Ch = linspace(-2.7, 5.0, length(S03(n).(['Chs_R' num2str(i)])));
    plot(time_Ch, S03(n).(['Chs_R' num2str(i)]), 'Color', [0.5 0.5 0.5]);
    hold on;
    plot(time_Ch, S03(n).(['Ch_R' num2str(i) 'avg']), 'k-', 'LineWidth', 2);
    ylim([Chs_R_min(i), Chs_R_max(i)]); % Set y-axis limits
    hold off;
end
%%
figure;
n = 4;
time = -2.7:1:5.0;


for i = 1:8
    subplot(4,2,i);
    time_Ch = linspace(-2, 5.7, length(S03(n).(['Chs_L' num2str(i)])));
    avg = S03(n).(['Ch_L' num2str(i) 'avg']);
    std_dev = S03(n).(['Ch_L' num2str(i) 'std']);
    
    
    % Plot average with thick black line
    %plot(time_Ch, avg,'Color', [148/255, 17/255, 0] ,'LineWidth', 3.5); 
    plot(time_Ch, avg, 'Color',[83/255, 27/255, 147/255] ,'LineWidth', 3.5);
    
    hold on;
    
    % Plot standard deviation as filled area
    
    % Plot standard deviation as filled area
    %fill([time_Ch, fliplr(time_Ch)], [avg + std_dev, fliplr(avg - std_dev)], [148/255, 17/255, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    fill([time_Ch, fliplr(time_Ch)], [avg + std_dev, fliplr(avg - std_dev)], [83/255, 27/255, 147/255], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    ylim([Chs_L_min(i), Chs_L_max(i)]); % Set y-axis limits
    if i == 3
       ylim([0, 1.5e-5]); % Set y-axis limits for subplot 3
    end
    xticks(-2:1:5.7);
    xticklabels({'','1','2','3','4','5','6','7'});
    hold off;
    
end

figure;
n = 3;
for i = 1:8
    subplot(4,2,i);
    time_Ch = linspace(-2, 5.7, length(S03(n).(['Chs_L' num2str(i)])));
    avg = S03(n).(['Ch_L' num2str(i) 'avg']);
    std_dev = S03(n).(['Ch_L' num2str(i) 'std']);
    
    % Plot average with thick black line
    %plot(time_Ch, avg, 'Color', [255/255, 147/255, 0], 'LineWidth', 3.5); 
    plot(time_Ch, avg, 'Color', [0, 84/255, 147/255], 'LineWidth', 3.5);
    hold on;
    
    %fill([time_Ch, fliplr(time_Ch)], [avg + std_dev, fliplr(avg - std_dev)], [255/255, 147/255, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    fill([time_Ch, fliplr(time_Ch)], [avg + std_dev, fliplr(avg - std_dev)], [0, 84/255, 147/255], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    ylim([Chs_L_min(i), Chs_L_max(i)]); % Set y-axis limits
    if i == 3
       ylim([0, 1.5e-5]); % Set y-axis limits for subplot 3
    end
    hold off;
    xticks(-2:1:5.7);
    xticklabels({'','1','2','3','4','5','6','7'});
    hold off;
end
%% Normalization
figure;
n = 1;
time = -2.7:1:5.0;


for i = 1:8
    subplot(4, 2, i);
    time_Ch = linspace(-2, 5.7, length(S03(n).(['Chs_L' num2str(i)])));
    
    % Get average and standard deviation
    avg = S03(n).(['Ch_L' num2str(i) 'avg']);
    std_dev = S03(n).(['Ch_L' num2str(i) 'std']);
    
    % Find the maximum value for normalization
    %max_val = max(max(S03(n).(['Chs_L' num2str(i)])));
    
    % Normalize the average and standard deviation
    avg_normalized = avg;
    std_dev_normalized = std_dev;
    
    % Plot average with thick black line
    plot(time_Ch, avg_normalized, 'Color', [83/255, 27/255, 147/255], 'LineWidth', 3.5);
    %plot(time_Ch, avg_normalized, 'Color', [255/255, 147/255, 0], 'LineWidth', 3.5);
    
    hold on;
    
    % Plot standard deviation as filled area
    fill([time_Ch, fliplr(time_Ch)], ...
         [avg_normalized + std_dev_normalized, fliplr(avg_normalized - std_dev_normalized)], ...
         [83/255, 27/255, 147/255], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    %fill([time_Ch, fliplr(time_Ch)], ...
         %[avg_normalized + std_dev_normalized, fliplr(avg_normalized - std_dev_normalized)], ...
         %[255/255, 147/255, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    % Set y-axis limits
    %ylim([Chs_L_min(i) / max_val, Chs_L_max(i) / max_val]);
    
    % Additional specific y-axis limit for subplot 3
    %if i == 3
       %ylim([0, 1.5e-5 / max_val]);
    %end
    
    % Set x-axis ticks and labels
    xticks(-2:1:5.7);
    xticklabels({'', '1', '2', '3', '4', '5', '6', '7'});
    ylim([0, 1]);
    
    hold off;
end
%%
figure;
n = 2;
time = -2.7:1:5.0;


for i = 1:8
    subplot(4, 2, i);
    time_Ch = linspace(-2, 5.7, length(S03(n).(['Chs_L' num2str(i)])));
    
    % Get average and standard deviation
    avg = S03(n).(['Ch_L' num2str(i) 'avg']);
    std_dev = S03(n).(['Ch_L' num2str(i) 'std']);
    
    % Find the maximum value for normalization
    %max_val = max(max(S03(n).(['Chs_L' num2str(i)])));
    
    % Normalize the average and standard deviation
    avg_normalized = avg;
    std_dev_normalized = std_dev;
    
    % Plot average with thick black line
    plot(time_Ch, avg_normalized, 'Color', [0, 84/255, 147/255], 'LineWidth', 3.5);
    %plot(time_Ch, avg_normalized, 'Color', [148/255, 17/255, 0], 'LineWidth', 3.5);
    
    hold on;
    
    % Plot standard deviation as filled area
    fill([time_Ch, fliplr(time_Ch)], ...
         [avg_normalized + std_dev_normalized, fliplr(avg_normalized - std_dev_normalized)], ...
         [0, 84/255, 147/255], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    %fill([time_Ch, fliplr(time_Ch)], ...
         %[avg_normalized + std_dev_normalized, fliplr(avg_normalized - std_dev_normalized)], ...
         %[148/255, 17/255, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    % Set y-axis limits
    %ylim([Chs_L_min(i) / max_val, Chs_L_max(i) / max_val]);
    
    % Additional specific y-axis limit for subplot 3
    %if i == 3
       %ylim([0, 1.5e-5 / max_val]);
    %end
    
    % Set x-axis ticks and labels
    xticks(-2:1:5.7);
    xticklabels({'', '1', '2', '3', '4', '5', '6', '7'});
    ylim([0, 1]);
    
    hold off;
end

%%


