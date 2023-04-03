%% generateRFPaths.m
% by Jim Sluijter
% This script checks the rigid foldability, and generates all rigidly
% foldable paths for all given geometries for a given RBM.
% The order of operation is created such that parallel computing can be
% utilized

function [X_EE,variables] = generateRFPaths(R2,validGeometries,Mode3,Mode4,inputAngles)

%allocate memory
temp_X_EE        = cell(size(validGeometries,2),1);         %allocate cells
temp_X_EE(:)     = {zeros(16,2)};                           %allocate cell content
temp_variables   = cell(size(validGeometries,2),1);
temp_variables(:)= {zeros(1,8)};

%
parfor i = 1:size(validGeometries,2)
    try
       thisVariables = [R2 validGeometries(:,i)' Mode3 Mode4];     %create vector with current variables [R2 a31 a32 a33 a41 a42 Mode3 Mode4]
       temp_X_EE{i} = calcPath(thisVariables,inputAngles);         %Calculate entire path for 
       temp_variables{i} = thisVariables;
    catch MExc
       temp_X_EE{i}      = []; %indicating error
       temp_variables{i} = []; %indicating error
       if ~(strcmp(MExc.message,'Not rigidly foldable!') || strcmp(MExc.message,'fold angle too large in vertex 3') || strcmp(MExc.message,'risk of singularity in vertex 4') || strcmp(MExc.message,'fold angle too large in vertex 4') )
           error(MExc.message)
       end
       continue
    end
end
nonemptyidx = ~cellfun(@isempty,temp_X_EE);
X_EE        = temp_X_EE(nonemptyidx); 
variables   = temp_variables(nonemptyidx);
