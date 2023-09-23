function [table]=final_ranked_predictions(connectivityMATRIX,genes,SIGN,fileNAME)

% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering 
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright. November 1, 2016.

if nargin < 2
    error('*** More input arguments needed. Please upload the data! ***');
end

if nargin < 3 || isempty(SIGN)
    SIGN= 0; % By default no sign info
end



% gene i=target; gene j=source;
numGENES=length(genes);

%Interactions
interactions=connectivityMATRIX(:);

% %Cutting
% idxZEROS=find(interactions~=0);
% interactions=interactions(idxZEROS);

%Edges
edges={};

if SIGN
    for i=1:length(interactions)
        if interactions(i)<0
            edges=[edges; 'repression'];
        else if interactions(i)>0
                edges=[edges; 'activation'];
            else
                edges=[edges; 'no regulation'];
            end
        end
    end
else
    for i=1:length(interactions)
        if interactions(i)==0
            edges=[edges; 'no regulation'];
        else 
            edges=[edges; 'activation/repression'];
        end
    end
end
% interactions=num2cell(abs(interactions));
interactions=num2cell(abs(interactions));
%Column of sources
sourceGENES=repmat(genes',1,numGENES)';
sourceGENES=sourceGENES(:);
% sourceGENES=sourceGENES((idxZEROS));
%Column of targets
targetGENES=repmat(genes,1,numGENES)';
targetGENES=targetGENES(:);
% targetGENES=targetGENES((idxZEROS));

titleCOLUMNS={'SourceGENES' 'TargetGENES' 'Interaction' 'Edges'};
excelTABLE=[sourceGENES targetGENES interactions edges];
table=cell2table(excelTABLE,'VariableNames',titleCOLUMNS);
table=sortrows(table, -3);
if nargin < 4 || isempty(fileNAME)
    fprintf('*** The results will not be saved! ***\n\n\n')
else
    if exist(fileNAME, 'file')
        delete(fileNAME);
    end
    writetable(table,fileNAME)
end
end