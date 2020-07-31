% if, whichone is not defined, it returns 0/1 (false/true) array 
% getDATA(data,strs,variable,value,whichone)

function[foundData] = getDATA(data,strs,variable,value,whichone)

tt = data(:,findCol(strs,variable{1}))==value{1};

if length(variable)>1
    for i=2:length(variable)
        tt = tt & data(:,findCol(strs,variable{i}))==value{i};
    end
end

if nargin==4
    foundData = tt;
else
    foundData = data(tt==1, whichone); 
end

return