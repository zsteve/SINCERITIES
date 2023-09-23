function [pred_lambda_min]=SINCERITIES_final(single_cell_data,time,numGENES,num_time_points,options,family,idx_lambda_min,noDIAG)

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
        [h(ti,gi),~,DISTANCE_matrix(ti,gi)]=kstest2(p1,p2);
    end
end

% figure;
% plot(DISTANCE_matrix)

%% Normalization
deltaT=time(2:end)-time(1:end-1);
DISTANCE_matrix_normed=DISTANCE_matrix./repmat(deltaT,1,size(DISTANCE_matrix,2));
DISTANCE_matrix=DISTANCE_matrix_normed;
X_matrix=DISTANCE_matrix(1:num_time_points-2,:);
%% Generate Y and X_matrix for glmnet

alphas=0;

%% GLMNET
pred_lambda_min=[];
for gi=1:numGENES
    if noDIAG==1
        options.exclude=gi;
    end
    
    options.alpha=alphas; %LASSO, for elasticnet set between 0-1
    Y_vector=DISTANCE_matrix(2:num_time_points-1,gi);
    CV_results_final = glmnet(X_matrix,Y_vector,family,options);
    pred_lambda_min(:,gi)=CV_results_final.beta(:,idx_lambda_min(gi));
end

end

