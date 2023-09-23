function [adj_matrix,DISTANCE_matrix_train]=SINCERITIES_PLUS(DATA,noDIAG,SIGN,CV_nfolds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             SINCERITIES_PLUS                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% SINCERITIES_PLUS is a novel computational method for inferring
% gene regulatory network (GRN) from time-stamped cross-sectional
% single-cell expression data.
% 

% SINCERITIES_PLUS function is an extension of the default version of
% SINCERITIES, designed for single cell datasets with fewer than five time
% points.
% DEFAULT regularization regression strategy: RIDGE
% DEFAULT distribution distance: KS (Kolmogorov-Smirnov) 
% 
% [adj_matrix,DISTANCE_matrix_train]=SINCERITIES_PLUS(DATA)
%
% DATA, a 1 by 1 structure containing the following information:
% - DATA.singleCELLdata: 1 by n cell array, where n is the number of
% capture time points. DATA.singleCELLdata{k} is a m by s_k matrix
% containing observed expression levels of m genes in s_k single cells.
% - DATA.totDATA: S by m matrix, where S is the total number of single
% cells (i.e., S=s_1+s_2+...+s_n where n the number of capture time points)
% and m is the number of genes.
% - DATA.time: n by 1 vector containing the cell capture time points or
% time-stamps).
%
% [adj_matrix,DISTANCE_matrix_train]=SINCERITIES_PLUS(DATA,noDIAG)
%
% noDIAG: this parameter selects whether the auto-regulatory edge is
% inferred
% 0- GRN contains no auto-regulatory edge (* DEFAULT *)
% 1- GRN contain auto-regulatory edge
%
% [adj_matrix,DISTANCE_matrix_train]=SINCERITIES_PLUS(DATA,noDIAG,SIGN)
% SIGN: this parameter selects whether the sign / mode of the gene
% regulations is inferred
% 0- for unsigned GRN
% 1- for signed GRN (* DEFAULT *)
% SINCERITIES uses partial correlation analysis where a positive (negative)
% correlation is taken as an indication of activation (repression).
%
% [adj_matrix,DISTANCE_matrix_train]=SINCERITIES_PLUS(DATA,noDIAG,SIGN,CV_n_folds)
% CV_n_folds defines a partition of the data into CV_n_folds disjoint
% subsets for the cross validation
% CV_n_folds=5; (* DEFAULT *)
%
% OUTPUTS:
%
% -adj_matrix: m by m matrix containing the weights of regulatory edges.
% The larger adj_matrix(i,j) indicates higher confidence that the
% corresponding edge exists (i.e., gene i regulating gene j).
% -DISTANCE_matrix: n-1 by m matrix containing the (normalized)
% distribution distance (DD) computed during the network inference.
%
% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright. Apr 20, 2017.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add dependecies
f=filesep;
[curr_path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(curr_path))

%% CHECK INPUT ARGUMENTS
if nargin < 1
    error_CV('*** More input arguments needed. Please upload the data! ***');
end

if nargin < 2 || isempty(noDIAG)
    noDIAG = 0; % Default inference considers self-regulatory edges
end

if nargin < 3 || isempty(SIGN)
    SIGN = 1;
end

if nargin < 4 || isempty(CV_nfolds)
    CV_nfolds = 5;
end

%% Initialization
single_cell_data=DATA.singleCELLdata;
time=DATA.time;
numGENES=size(single_cell_data{1},1);
num_time_points=length(time);

%% Catch error
if num_time_points<3
    error('** The data must contain at least 3 time points **')
end

%% GLMNET options
options=glmnetSet;
options.cl=[0;Inf];  %ONLY POS COEFFICIENT ESTIMATION
options.standardize=true;
%     options.intr=0;
options.lambda = logspace(-2,2,100);
% options.intr=0;
family=[];

%% K-fold CV
K=CV_nfolds;
for cv=1:size(single_cell_data,2)
    N=size(single_cell_data{cv},2);
    indices{cv} = crossvalind('Kfold', N, K);
end
for cross=1:K
    for cv=1:size(single_cell_data,2)
        test = (indices{cv} == cross); train = ~test;
        data4test{cv}=single_cell_data{cv}(:,test);
        data4train{cv}=single_cell_data{cv}(:,train);
    end
    
    %% ********** DISTRIBUTION DISTANCE train **************
    
    DISTANCE_matrix_train=zeros(num_time_points-1,numGENES);
    h=zeros(num_time_points-1,numGENES);
    totalDATA=data4train{1}';
    
    for ti=1:num_time_points-1
        totalDATA=[totalDATA; data4train{ti+1}'];
        data_ti=data4train{ti};
        data_ti_plus1=data4train{ti+1};
        for gi=1:numGENES %parfor gi=1:numGENES
            p1=data_ti(gi,:);
            p2=data_ti_plus1(gi,:);
            [h(ti,gi),~,DISTANCE_matrix_train(ti,gi)]=kstest2(p1,p2);
        end
    end
    
    % Normalization
    deltaT=time(2:end)-time(1:end-1);
    DISTANCE_matrix_train=DISTANCE_matrix_train./repmat(deltaT,1,size(DISTANCE_matrix_train,2));
    X_matrix=DISTANCE_matrix_train(1:num_time_points-2,:);
    %% ********** DISTRIBUTION DISTANCE test **************
    
    DISTANCE_matrix_test=zeros(num_time_points-1,numGENES);
    h=zeros(num_time_points-1,numGENES);
    totalDATA=data4test{1}';
    
    for ti=1:num_time_points-1
        totalDATA=[totalDATA; data4test{ti+1}'];
        data_ti=data4test{ti};
        data_ti_plus1=data4test{ti+1};
        for gi=1:numGENES %parfor gi=1:numGENES
            p1=data_ti(gi,:);
            p2=data_ti_plus1(gi,:);
            [h(ti,gi),~,DISTANCE_matrix_test(ti,gi)]=kstest2(p1,p2);
        end
    end
    
    % Normalization
    deltaT=time(2:end)-time(1:end-1);
    DISTANCE_matrix_test=DISTANCE_matrix_test./repmat(deltaT,1,size(DISTANCE_matrix_test,2));
    X_matrix_test=DISTANCE_matrix_test(1:num_time_points-2,:);
    
    %% Generate Y and X_matrix for glmnet
    alphas=0; % RIDGE
    
    % GLMNET
    pred_lambda_min=[];
    for gi=1:numGENES
        if noDIAG==1
            options.exclude=gi;
        end
        
        options.alpha=alphas; %LASSO, for elasticnet set between 0-1
        Y_vector=DISTANCE_matrix_train(2:num_time_points-1,gi);
        CV_results = glmnet(X_matrix,Y_vector,family,options);
        Y_vector_test=DISTANCE_matrix_test(2:num_time_points-1,gi);
        
        for lambdacount=1:length(CV_results.lambda)
            beta_lambda=CV_results.beta(:,lambdacount);
            error_CV(cross,gi,lambdacount)=sum((Y_vector_test-X_matrix_test*beta_lambda).^2);
        end
    end
    
end

mean_error_CV=squeeze(mean(error_CV,1));
standard_error_mean=squeeze(std(error_CV,1))/sqrt(K);

% Lambda_min
[min_mean_error_CV,idx_lambda_min]=min(mean_error_CV,[],2);
lambda_min=CV_results.lambda(idx_lambda_min);
% Lambda 1SE
for gi=1:numGENES
    min_plu_1SE=mean_error_CV(gi,idx_lambda_min(gi))+standard_error_mean(gi,idx_lambda_min(gi));
    idx_lambda_1SE(gi) = find(mean_error_CV(gi,:)<=min_plu_1SE,1);
end
lambda_1SE=CV_results.lambda(idx_lambda_1SE);

% figure
% plot(-log10(CV_results.lambda),mean_error_CV')
% xlabel('-Log(lambda)')
% ylabel('Mean-squared error')
% grid on
% hold on
% for gi=1:numGENES
%     plot(-log10(lambda_min(gi)),mean_error_CV(gi,idx_lambda_min(gi)),'r*');
%     plot(-log10(lambda_1SE(gi)),mean_error_CV(gi,idx_lambda_1SE(gi)),'b*');
% end

pred_lambda=SINCERITIES_final(single_cell_data,time,numGENES,num_time_points,options,family,idx_lambda_min,noDIAG);


if SIGN
    %% Partial correlation analysis
    
    [parcorr_matrix]=partialcorr(DATA.totDATA,'type','Spearman');
    pred_lambda=pred_lambda.*(sign(parcorr_matrix));
    
    
end

adj_matrix=pred_lambda;

end