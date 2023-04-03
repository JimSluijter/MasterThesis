%% EvaluatePaths2
% by Jim Sluijter This script performs further evaluations on the resulting
% paths, calling functions to discard all geometries with sector angles
% under 10 degrees, and all geometries where the dihedral angles form a
% risk.
% This script has been created after the generation of the paths, and is
% therefore seperate from the other evaluatePaths.m script.

% runID = '230113_05';
runID = '230201_10';
load(sprintf('Results/%s/%s_D_analyzed.mat',runID,runID),'feasible');
input = generateInputProfile([15 80],16);
%alocate memory for new feasible structure
% new_feasible = zeros(size(feasible));
k = 1;

for i = 1:length(feasible)
    variables = feasible{i}.variables;
    
    %1. Filter out all geometries with sector angles below 10 degrees
    if any(variables(2:6) < 10)
        continue
    end
    
    %2. Filter out all rows with high fold angles
    dihAngleStatus = filterDihedralAngles10(variables,input);
    if dihAngleStatus ~= 0
        continue
    end
    
    %save survivors
    new_feasible{k} = feasible{i};
    k = k+1;
end     
% delete empty cells in new_feasible and rename to feasible

clear feasible
feasible = new_feasible;


save(sprintf('Results/%s/%s_E_extra_requirements.mat',runID,runID),'feasible');