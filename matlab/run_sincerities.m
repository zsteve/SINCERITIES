clear all
close all
clc

f=filesep;
[curr_path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(curr_path))

%% *** Data loading ***
kounoFILEname=fullfile(pwd,"X.xlsx");
DATA=uploading(kounoFILEname);

%% *** SINCERITIES ***

% Parameter settings:
distance=1; %Kolmogorov-Smirnov
method=1; %Ridge regression
noDIAG=1; %GRN inference includes autoregulatory edges
SIGN=0; %Predict unsigned GRN
[adj_matrix,DISTANCE_matrix]=SINCERITIES(DATA,distance,method,noDIAG,SIGN);

%%
adj_matrix_norm=adj_matrix/max(max(adj_matrix)); % Normalization
filename=fullfile(pwd,'A');
[table]=final_ranked_predictions(adj_matrix_norm,DATA.genes,SIGN,filename);
exit
