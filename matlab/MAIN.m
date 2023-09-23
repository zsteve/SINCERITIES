%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      SINCERITIES MAIN SCRIPT                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ** Please prepare the data in an Excel sheet using the format below **
% Data: s-by-m+1 matrix, where s is the total number of observations/single 
% cells and m is the number of genes. The first m columns contain 
% the expression level of each m genes, and the last column contains the
% the time-stamps.
%
% Two data formats are accepted:
% 
% A) with row header
% ------------------------------------
% Gene1  Gene2  Gene3  ... Genej  Time
%  27     80     56    ...  69      0
%  73     20     90    ...  45      0
%   .     .      .     ...  .       .
%   .     .      .     ...  .       .
%   .     .      .     ...  .       .
% ------------------------------------
% 
% B) without row header
% ------------------------------------
%  27     80     56    ...  69      0
%  73     20     90    ...  45      0
%   .     .      .     ...  .       .
%   .     .      .     ...  .       .
%   .     .      .     ...  .       .
% ------------------------------------
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

f=filesep;
[curr_path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(curr_path))

%% *** Data loading ***

DATA=uploading();

%% *** SINCERITIES ***

% Parameter settings:
% distance: this parameter selects the distribution distance
% 1- for KS (Kolmogorov-Smirnov)  (* DEFAULT *)
% 2- for CM (Cramer-von Mises)
% 3- for AD (Anderson-Darling)
%
% method: this parameter selects the regularization regression strategy
% 1- for RIDGE  (* DEFAULT *)
% 2- for ELASTIC-NET with automatic detection of optimal alpha parameter
% 3- for LASSO
% 4- for ELASTIC-NET with manual selection of alpha parameter
%
% noDIAG: this parameter selects whether the auto-regulatory edge is
% inferred
% 0- GRN contains auto-regulatory edges (* DEFAULT *)
% 1- GRN contains no auto-regulatory edge
%
% SIGN: this parameter selects whether the sign / mode of the gene
% regulations is inferred
% 0- for unsigned GRN
% 1- for signed GRN (* DEFAULT *)

distance=1;
method=1;
noDIAG=0;
SIGN=1;

[adj_matrix,DISTANCE_matrix]=SINCERITIES(DATA,distance,method,noDIAG,SIGN);

% %% For data with less than 5 time points please use:
% CV_n_folds=10;
% [adj_matrix,DISTANCE_matrix]=SINCERITIES_CROSS_VALIDATION(DATA,noDIAG,SIGN,CV_n_folds);

%% Final ranked list of regulatory edges
adj_matrix=adj_matrix/max(max(adj_matrix));
filename=fullfile(pwd,'Results','GRNprediction');

%% Final ranked list of predicted edges
[table]=final_ranked_predictions(adj_matrix,DATA.genes,filename);
disp(table);