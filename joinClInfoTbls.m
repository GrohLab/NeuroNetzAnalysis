function [clInfoOut] = joinClInfoTbls(clInfoA, clInfoB)
%JOINCLINFOTBLS tries to homogenize the variables in each table such that
%they can be merged. This function is not designed to merge any kind of
%tables but only those which Phy outputs.
%   Basically it places table A on top of table B.
%   [clInfoOut] = joinClInfoTbls(clInfoA, clInfoB);
% Emilio Isaias-Camacho @ GrohLab 2020

% Phy different names for the same variable
phyNames = ["cluster_id","id";...
    "ch","channel";...
    "sh", "shank";...
    "fr", "firing_rate"];

popVarNames = clInfoA.Properties.VariableNames;
curVarNames = clInfoB.Properties.VariableNames;
Ncv = numel(curVarNames); % Npv = numel(popVarNames); 
curInPopFlags = ismember(popVarNames, curVarNames);
% Modifying the current cluster info to match the population table
missingInCurSub = ~curInPopFlags;
missingNames = popVarNames(missingInCurSub);
% Are the names which are missing a user addition or a phy update?
for mn = missingNames
    addnv = true;
    phyFlag = ismember(phyNames, mn);
    if any(phyFlag(:)) % Is this name a phy update?
        addnv = false;
        varSub = (1:size(phyFlag,1))*any(phyFlag, 2); % Which variable
        phyVer = any(phyFlag, 1)*[1;2];
        fprintf(1, 'Substitution of %s for %s\n',...
            phyNames(varSub, 3-phyVer), phyNames(varSub, phyVer));
        posInCur = strcmp(curVarNames, phyNames(varSub, 3-phyVer));
        if any(posInCur)
            clInfoB.Properties.VariableNames(posInCur) =...
                cellstr(phyNames(varSub, phyVer));
        else
            fprintf(1, 'No substitution, rather ')
            addnv = true;
        end
    end
    if addnv% or a user added variable
        fprintf(1, 'Adding %s to the current table\n', mn{1});
        tempTable = groupsummary(clInfoA, mn);
        [~, mstSub] = sortrows(tempTable, 'GroupCount', 'descend');
        fillVal = tempTable{mstSub(1),mn};
        clInfoB = addvars(clInfoB, repmat(fillVal, size(clInfoB,1), 1),...
            'NewVariableNames', mn);
    end
end
curVarNames = clInfoB.Properties.VariableNames;
[~, popInCurSubs] = ismember(popVarNames, curVarNames);
% Modifying the population cluster info to match the current table
missingInPopSubs = setdiff(1:Ncv, popInCurSubs);
missingNames = curVarNames(missingInPopSubs);
% I'm assuming that the missing variables are only user added
% names.
for mn = missingNames
    fprintf(1, 'Adding %s to the population table\n', mn{1});
    tempTable = groupsummary(clInfoB, mn);
    [~, mstSub] = sortrows(tempTable, 'GroupCount', 'descend');
    fillVal = tempTable{mstSub(1),mn};
    clInfoA = addvars(clInfoA, repmat(fillVal, size(clInfoA,1),1),...
        'NewVariableNames', mn);
end
popVarNames = clInfoA.Properties.VariableNames;
[curInPopFlags, popInCurSubs] = ismember(popVarNames, curVarNames);
if ~all(curInPopFlags) || ~isempty(setdiff(1:Ncv, popInCurSubs))
    fprintf(1, 'Predicting an error!\n')
end
clInfoOut = cat(1, clInfoA, clInfoB);
end

