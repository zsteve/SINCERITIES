%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 EXAMPLE 2: THP-1 single cell data                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

f=filesep;
[curr_path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(curr_path))
% addpath(fullfile(pwd,f,'THP1 data'))

%% *** Data loading ***
kounoFILEname=fullfile(curr_path,'THP1 data','THP1_single_cell_data_EXCEL_no6_24_72_96.xlsx');
DATA=uploading(kounoFILEname);

% %% ************** SUBNET SELECTION BEFORE NET INFERENCE *************
% 
% [~,subGENES]= xlsread(fullfile(pwd,'THP1 data','SUBNET2_tomaru.xlsx'));
% [~,idxSUBgenes_first]=ismember(subGENES,DATA.genes);
% DATA.numGENES=length(idxSUBgenes_first);
% DATA.genes=DATA.genes(idxSUBgenes_first);
% DATA.totDATA=DATA.totDATA(:,idxSUBgenes_first);
% for i=1:DATA.num_time_points
% DATA.singleCELLdata{i}=DATA.singleCELLdata{i}(idxSUBgenes_first,:);
% end

%% *** SINCERITIES ***

% Parameter settings:
noDIAG=0; %GRN inference includes autoregulatory edges
SIGN=0; %Predict unsigned GRN
CV_nfolds=10;
[adj_matrix,DISTANCE_matrix]=SINCERITIES_PLUS(DATA,noDIAG,SIGN,CV_nfolds);

%% *** SUBNETWORK EXTRACTION FOR OVERLAPPING TFs ***
[~,subGENES]= xlsread(fullfile(curr_path,'THP1 data','SUBNET2_tomaru.xlsx'));
[~,idxSUBgenes]=ismember(subGENES,DATA.genes);
idxSUBgenes(idxSUBgenes==0)=[];
DATA.numGENES=length(idxSUBgenes);
DATA.genes=DATA.genes(idxSUBgenes);
adj_matrix=adj_matrix(idxSUBgenes,idxSUBgenes);
for i=1:DATA.num_time_points
    DATA.singleCELLdata{i}=DATA.singleCELLdata{i}(idxSUBgenes,:);
end
DATA.totDATA=DATA.totDATA(:,idxSUBgenes);

%% *** REFERENCE GRN from RNAi EXPERIMENTS (from Tomaru et al.) ***

[type_regulation,netINFO]= xlsread(fullfile(curr_path,'THP1 data','tomaru2.xlsx'));
adj_ref=zeros(DATA.numGENES);
for i=1: size(netINFO,1)
    [exists, idxGENEsource]=ismember(netINFO(i,1), DATA.genes);    
    if exists==1        
        [~, idxGENEtarget]=ismember(netINFO(i,3:end), DATA.genes);
        idxGENEtarget(idxGENEtarget==0)=[];
        if SIGN==1 %Assign edge signs to the reference GRN if SIGN = 1 
            adj_ref(idxGENEsource,idxGENEtarget) = type_regulation(i);
        else
            adj_ref(idxGENEsource,idxGENEtarget) = 1;
        end
    end
end
netINFO(:,2)=[];
tomaruGENES=unique(reshape(netINFO,size(netINFO,1)*size(netINFO,2),1)); 
[found,iii]=ismember('ND',tomaruGENES); 
if found==1
    tomaruGENES(iii)=[];
end

%% Final ranked list
adj_matrix_norm=adj_matrix/max(max(adj_matrix)); % Normalization
filename=fullfile(curr_path,'Results','prediction4THP1');
[table]=final_ranked_predictions(adj_matrix_norm,DATA.genes,SIGN,filename);
fprintf('*** Inferred regulatory relationships for THP-1 single cell data *** \n\n\n');
disp(table);

%% AUROC (x=fpr/1-specifity; y=recall/sensitivity) and AUPR (x=recall y=precision)
% Auto-regulatory edges are removed for AUROC and AUPR evaluation since
% RNAi experiments would not allow the identification of such edges. 
adj_matrix_norm=remove_diagonal(adj_matrix_norm); %Removing auto-regulatory edges
[AUC(1),AUC(2)]=auc_from_ranks_TC_sign(adj_matrix_norm,adj_ref,1000);
AUC