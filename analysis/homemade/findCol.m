function [whichone] = findCol(strs,variable)

if ~iscell(variable)
    whichone = find(ismember(strs,variable));
else
    whichone = [];
    for i=1:numel(variable)
        t_whichone = find(ismember(strs,variable{i}));
        whichone = [whichone t_whichone];
    end
end
return