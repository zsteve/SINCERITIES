function [adj_matrix,DISTANCE_matrix]=SINCERITIES(DATA,distance,method,noDIAG,SIGN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             SINCERITIES                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% SINCERITIES is a novel computational method for inferring
% gene regulatory network (GRN) from time-stamped cross-sectional
% single-cell expression data.

% [adj_matrix,DISTANCE_matrix]=SINCERITIES(DATA) 
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
% [adj_matrix,DISTANCE_matrix]=SINCERITIES(DATA,distance) 
%
% distance: this parameter selects the distribution distance
% 1- for KS (Kolmogorov-Smirnov)  (* DEFAULT *)
% 2- for CM (Cramer-von Mises)
% 3- for AD (Anderson-Darling)
% 4- for Mean expression difference
%
% [adj_matrix,DISTANCE_matrix]=SINCERITIES(DATA,distance,method) 
%
% method: this parameter selects the regularization regression strategy
% 1- for RIDGE  (* DEFAULT *)
% 2- for ELASTIC-NET with automatic detection of optimal alpha parameter
% 3- for LASSO
% 4- for ELASTIC-NET with manual selection of alpha parameter
%
% [adj_matrix,DISTANCE_matrix]=SINCERITIES(DATA,distance,method,noDIAG)
%
% noDIAG: this parameter selects whether the auto-regulatory edge is
% inferred
% 0- GRN contains no auto-regulatory edge (* DEFAULT *)
% 1- GRN contain auto-regulatory edge
%
%
% [adj_matrix,DISTANCE_matrix]=SINCERITIES(DATA,distance,method,noDIAG,SIGN)
% SIGN: this parameter selects whether the sign / mode of the gene
% regulations is inferred
% 0- for unsigned GRN
% 1- for signed GRN (* DEFAULT *)
% SINCERITIES uses partial correlation analysis where a positive (negative) 
% correlation is taken as an indication of activation (repression). 
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
% Copyright. November 1, 2016.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add dependecies
f=filesep;
[curr_path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(curr_path,f,'cmtest'))
addpath (fullfile(curr_path,'glmnet_matlab'))
addpath (fullfile(curr_path,'AnDarksamtest'))

%% CHECK INPUT ARGUMENTS
if nargin < 1
    error('*** More input arguments needed. Please upload the data! ***');
end

if nargin < 2 || isempty(distance)
    distance = 1; % Default distribution distance: KS
end

if nargin < 3 || isempty(method)
    method = 1; % Default regularization method: RIDGE
end

if nargin < 4 || isempty(noDIAG)
    noDIAG = 0; % Default inference considers self-regulatory edges
end

if nargin < 5 || isempty(SIGN)
    SIGN = 1; % Default inference considers self-regulatory edges
end

%% Initialization
single_cell_data=DATA.singleCELLdata;
time=DATA.time;
numGENES=size(single_cell_data{1},1);
num_time_points=length(time);

%% Catch error
if num_time_points<5
    error('** DATA with a number of time points < 5. Please run SINCERITIES_CROSS_VALIDATION function **')
end

%% GLMNET options
options=glmnetSet;
options.cl=[0;Inf];  %ONLY POS COEFFICIENT ESTIMATION
options.standardize=true;
%     options.intr=0;
%     option.lambda = 10.^linspace(-0.2,-6,100);
family=[];%'poisson';%[];  % Default from gaussian distribution
type=[];
parallel=[];

%% ********** DISTRIBUTION DISTANCE **************

DISTANCE_matrix=zeros(num_time_points-1,numGENES);
h=zeros(num_time_points-1,numGENES);
totalDATA=single_cell_data{1}';

for ti=1:num_time_points-1
    totalDATA=[totalDATA; single_cell_data{ti+1}'];
    data_ti=single_cell_data{ti};
    data_ti_plus1=single_cell_data{ti+1};
    for gi=1:numGENES %parfor gi=1:numGENES
        p1=data_ti(gi,:);
        p2=data_ti_plus1(gi,:);
        %         p1(p1==0)=[];
        %         p2(p2==0)=[];
        switch distance
            case 1
                
                [h(ti,gi),~,DISTANCE_matrix(ti,gi)]=kstest2(p1,p2);
                
            case 2
                
                [h(ti,gi),~,DISTANCE_matrix(ti,gi)]=cmtest2(p1,p2);
                
            case 3
                
                X=[p1';p2'];
                samples=[ones(length(p1),1); 2*ones(length(p2),1)];
                X=[X samples];
                [h(ti,gi),DISTANCE_matrix(ti,gi)]=AnDarksamtest(X);
            case 4
                [DISTANCE_matrix(ti,gi)]=abs(mean(p2)-mean(p1));     
            otherwise
                
                disp('**** ERROR: CHOOSE 1 for KS, 2 for CM, and 3 for AD! ****')
        end
        
        
    end
end

% figure;
% plot(DISTANCE_matrix)

%% Normalization
deltaT=time(2:end)-time(1:end-1);
DISTANCE_matrix_normed=DISTANCE_matrix./repmat(deltaT,1,size(DISTANCE_matrix,2));

%% Generate Y and X_matrix for glmnet

switch method
    case 1
        alphas=0;
    case 2
        alphas=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
    case 3
        alphas=1;
    case 4
        alphas=input(' *** Please input manually the alpha value (between 0 and 1): ');
    otherwise
        disp('**** ERROR: CHOOSE 1 for RIDGE, 2 for ELASTIC NET, 3 for LASSO! ****')
end
DISTANCE_matrix=DISTANCE_matrix_normed;
X_matrix=DISTANCE_matrix(1:num_time_points-2,:);

%% LOOCV settings
nfolds=size(X_matrix,1);
foldid=1:nfolds;
keep=true;
pred_lambda_min=[];
% pred_lambda_1se=[];

for gi=1:numGENES
    if noDIAG==1
        options.exclude=gi;
    end
    for test=1:length(alphas)
        options.alpha=alphas(test); %LASSO, for elasticnet set between 0-1
        Y_vector=DISTANCE_matrix(2:num_time_points-1,gi);
        CV_results = cvglmnet(X_matrix,Y_vector,family,options,type,nfolds,foldid,parallel,keep,0);
        elasticNETresults{test}=CV_results;
        cvERROR(test)=min(CV_results.cvm);
    end
    %             figure
    %             plot(alphas,cvERROR)
    [~,idxMINerror]=min(cvERROR);
    bestALPHA=alphas(idxMINerror);
    alphaCHOSEN(gi)=bestALPHA;
    i_lambda_min=find(elasticNETresults{idxMINerror}.lambda==elasticNETresults{idxMINerror}.lambda_min);
    %     i_lambda_1se=find(elasticNETresults{idxMINerror}.lambda==elasticNETresults{idxMINerror}.lambda_1se);
    beta_lambda_min=elasticNETresults{idxMINerror}.glmnet_fit.beta(:,i_lambda_min);
    %     beta_lambda_1se=elasticNETresults{idxMINerror}.glmnet_fit.beta(:,i_lambda_1se);
    pred_lambda_min(:,gi)=beta_lambda_min;
    %     pred_lambda_1se(:,gi)=beta_lambda_1se;
    
end
pred_lambda=pred_lambda_min;

if SIGN
    %% Partial correlation analysis
    
    %     for i=1:num_time_points
    %         data_3D(:,:,i)=single_cell_data{i};
    %     end
    %     data_time_series=permute(data_3D,[2 3 1]);
    %
    [parcorr_matrix]=partialcorr(DATA.totDATA,'type','Spearman');
    
    pred_lambda=pred_lambda.*(sign(parcorr_matrix));
    
    
    %     % CORRELATION for each time point
    %
    %     for i=1:num_time_points
    %         [parcorr_matrix(:,:,i),~]=partialcorr(squeeze(data_3D(:,:,i))','type','Spearman');
    %
    %     end
    %
    %     [~,maxPARidx]=max(abs(parcorr_matrix),[],3);
    %     for ii=1:numGENES
    %         for jj=1:numGENES
    %             maxPARcorr(ii,jj)=parcorr_matrix(ii,jj,maxPARidx(ii,jj));
    %         end
    %     end
    %
    %     pred_lambda=pred_lambda.*(sign(maxPARcorr));
end

adj_matrix=pred_lambda;

end