%% HEATMAPS
%by Jim Sluijter
% used to create the heatmaps for all different criteria
clear
runID = '230113_06'; rbms = 1:4;
% runID = '230201_08'; rbms = 3;
load(sprintf('Results/%s/%s_D_analyzed.mat',runID,runID),'feasible','tbl','splittbl');


%%
% Number of feasible paths
% f(1,:) = 
makeHeatmaps(splittbl,rbms,"NumberOfFeasiblePaths",runID);
% Stride length
makeHeatmaps(splittbl,rbms,"strideLength",runID);
% Step height
makeHeatmaps(splittbl,rbms,"stepHeight",runID);
% Ground contact ratio
makeHeatmaps(splittbl,rbms,"contactRatio",runID);
% Error_height
makeHeatmaps(splittbl,rbms,"error_Z",runID);
% Error_speed
makeHeatmaps(splittbl,rbms,"error_speed",runID);
