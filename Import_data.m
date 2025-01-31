
clear; clc;
%%
params.exo = {'noexo', 'exo'}; % Never change the order
params.leg = {'right'}; % Change based on the leg needed to be analyzed
params.step = {'normal', 'max'}; % Never change the order
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
params.MotVars = {'RHipAngle','RKneeAngle','RAnkleAngle', ... % right to left, angles,
       'LHipAngle','LKneeAngle','LAnkleAngle', ...
       'ThoraxAngle','R1Met','L1Met',  ...
       'RHipMoment','RKneeMoment','RAnkleMoment', ... % Moments
       'LHipMoment','LKneeMoment','LAnkleMoment', ...
       'RGRF', 'LGRF','L_ASIS','R_ASIS', 'RLoadCell', 'LLoadCell'};     % External forces
datafolder = 'C:\Users\mpasav2\OneDrive - University of Illinois Chicago\Desktop\All_particpants_GI_study\All Subjects_Data\S05';
params.Motion_datafolder = [datafolder '\QTM_mat\']; % folder contatining Qualysis trial data
params.EMG_datafolder = [datafolder '\EMG_mat\']; % folder containing EMG trial data
params.Onset_datafolder = [datafolder '\onsetts\'];  % folder containing Onset data
params.trialnum = {
1:3, ... % Baseline trials
4:10, ... % WO R Norm
11:17, ... % W0 R Long
18:23, ... % W R Norm
24:30, ... % W R Long
};
% trial numbers for each condition
params.filenum = {1:3, ... % Baseline trials
  21:27, ... % WO R Norm
  8:14, ... % WO R Long
  15:20, ... % W R Norm
  1:7, ... % W R Long
};
params.musnames = {'Gmax','Gmed','Adl','RF','VL','BF','TA','MG'};



%% Importing 
S01 = EMG_import(params);
S02 = Mocap_import(params,S01);
S03 = onset_import(params,S02);
%S03(end) = [];

%% Filtering data
var = 'LHipAngle(1,:)'; % choose the variable
%var = 'LLoadCell';
baseline_correction = 1;    % choose whether to perform baseline correction or not 

params.filenum = {1:3, ... % Baseline trials
  21:27, ... % WO R Norm
  8:14, ... % WO R Long
  15:20, ... % W R Norm
  1:7, ... % W R Long
};
params.musnames = {'Gmax','Gmed','Adl','RF','VL','BF','TA','MG'};

processed_data = preprocessing(params, var, S03, baseline_correction);
%save('C:\Users\mpasav2\OneDrive - University of Illinois Chicago\Desktop\All Subjects_Data/S01_allTrials.mat', 'S03');
%save('C:\Users\mpasav2\OneDrive - University of Illinois Chicago\Desktop\All Subjects_Data/S01_Processed_LhipTorquesY.mat', 'processed_data');

%% Plotting data
opt = 'avg'; % choose whether to plot individual trials or avg+/-std; avg or ind
plotting(params, processed_data, var, opt);

%% Appending torques into the dataset
S04 = torque_compute(params, S03);
save('C:\Users\mpasav2\OneDrive - University of Illinois Chicago\Desktop\All_particpants_GI_study/All Subjects_Data/C01_allTrials_corrected.mat', 'S04');

%% saving data for summary plot
%{
Subject = 'S01';
filename_frontal_exo = [Subject 'baseline_Y_frontal_max_normal_exo_' params.leg{1}];
filename_frontal_noexo = [Subject 'baseline_Y_frontal_max_normal_noexo_' params.leg{1}];

filename_sagittal_exo = [Subject 'baseline_Y_sagittal_max_normal_exo_' params.leg{1}];
filename_sagittal_noexo = [Subject 'baseline_Y_sagittal_max_normal_noexo_' params.leg{1}];
if var(end-3) == '2'
    hippos = processed_data(1).raw_data; 
    hipos = processed_data(2).raw_data; 
    save(filename_frontal_noexo,'hippos','hipos');

    clearvars hipos hippos
    hippos = processed_data(3).raw_data; 
    hipos = processed_data(4).raw_data; 
    save(filename_frontal_exo,'hippos','hipos');
elseif var(end-3) == '1'
    hippos = processed_data(1).raw_data; 
    hipos = processed_data(2).raw_data; 
    save(filename_sagittal_noexo,'hippos','hipos');

    clearvars hipos hippos
    hippos = processed_data(3).raw_data; 
    hipos = processed_data(4).raw_data; 
    save(filename_sagittal_exo,'hippos','hipos');

end
%}
