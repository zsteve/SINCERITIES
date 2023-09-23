function [DATA]=uploading(filename)

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
% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering 
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright. November 1, 2016.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n**** Please select excel file containing data ****\n\n')

if nargin < 1 || isempty(filename)
    [FileName,PathName,FilterIndex] = uigetfile('*.*');
    filename=strcat(PathName, FileName);
end

[NUM,TXT]= xlsread(filename);
NUM(isnan(NUM(:,1)),:)=[]; % REMOVE ROW WITH AT LEAST ONE NaN
totDATA=NUM(:,1:end-1);
timeline=NUM(:,end);
DATA.time=unique(timeline);
DATA.num_time_points=length(DATA.time);
sortTOTdata=[];
sortTIMELINE=[];
for k=1:DATA.num_time_points
    I=find(timeline==DATA.time(k));
    cutDIMENSION(k)=length(I);
    sortTOTdata=[sortTOTdata; totDATA(I,:)];
    sortTIMELINE=[sortTIMELINE; timeline(I)];
end
DATA.totDATA=sortTOTdata;
DATA.timeline=sortTIMELINE;
%%%%%%%%
DATA.totDATA(isnan(totDATA)) = 0 ;
%%%%%%%%
DATA.numGENES=size(DATA.totDATA,2);
if isempty(TXT)
    for i=1:DATA.numGENES
        DATA.genes{i}=sprintf( '%s %i ', 'Gene ', i);
    end
else
    DATA.genes=TXT(1,1:end-1);
end
data= mat2cell(DATA.totDATA',DATA.numGENES,cutDIMENSION);% now rows=genes and columns=cells
DATA.singleCELLdata=data;
end