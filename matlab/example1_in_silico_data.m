%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               EXAMPLE 1: in silico single cell data                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

f=filesep;
[curr_path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(curr_path))

%% *** Data loading ***

load('20_nets_10genes_8UNEVENtime_sigma01B_no_initial_points2.mat')
DATA.time=time_points';
DATA.numGENES=n;
for numEXAMPLES=1:size(networks,3)
    % Data Preparation
    data_time_series=data_tot_array{numEXAMPLES};
    for i=1:num_time_points
        DATA.singleCELLdata{i}=squeeze(data_time_series(:,i,:))';
    end
    for i=1:DATA.numGENES
        DATA.genes{i}=sprintf( '%s %i ', 'Gene ', i);
    end
    DATA.totDATA=[];
    for i=1:num_time_points
        DATA.totDATA=[DATA.totDATA; squeeze(data_time_series(:,i,:))];
    end
    
    %% *** SINCERITIES ***
    
    % Parameter settings:
    distance=1; %Kolmogorov-Smirnov distance
    method=1; %Ridge regression
    noDIAG=1; %Assume GRN contains no autoregulatory edges 
    SIGN=1; %Predict signed GRN
    
    [adj_matrix,DISTANCE_matrix]=SINCERITIES(DATA,distance,method,noDIAG,SIGN);
    
    %% *** Performance Evaluation ***
    
    % Gold standard GRN 
    a=squeeze(networks(:,:,numEXAMPLES)); %GRNs from GeneNetWeaver
    a=remove_diagonal(a); %Remove auto-regulations as in silico data were generated without such edges
    if SIGN==0
        a(a~=0)=1;  %Ignoring edge signs if SIGN option was set to 0
    end
    
    %% Final ranked list
    % Normalization
    adj_matrix_norm=adj_matrix/max(max(adj_matrix));
    
    %% AUROC (x=fpr/1-specifity; y=recall/sensitivity) and AUPR (x=recall y=precision)
    [AUROC(numEXAMPLES),AUPR(numEXAMPLES)]=auc_from_ranks_TC_sign(adj_matrix_norm,a,1000);
    txtname=sprintf( '%s%i ', 'results4insilicoNETWORK', numEXAMPLES);
    filename=fullfile(pwd,'Results',txtname);
    [table]=final_ranked_predictions(adj_matrix_norm,DATA.genes,SIGN,filename);
    %     disp(table);
end
AUC=[AUROC' AUPR'];
mean_std=[mean(AUC(:,1)) mean(AUC(:,2)); std(AUC(:,1)) std(AUC(:,2))]
AUC=[AUC; mean_std];

%% Example: inferred regulatory relationships for network 20
fprintf('    *** Inferred regulatory relationships for network 20 *** \n\n\n');
disp(table);


