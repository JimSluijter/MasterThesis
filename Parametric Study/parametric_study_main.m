%% Paramteric sweep %% 
% by Jim Sluijter

% This script is used to call the parametric study.
% The generation of valid geometries, the generation of rigidly foldable
% paths, and their evaluation are called seperately for small portions of
% the entire data set to ensure good memory management. 
%%
clear
%% RUN ID
% used to identify the stored results
metaData.runID   = '230201_08_TT';

%%Current calculation: 
R2Range          = 1:7;
metaData.R2Range = R2Range;
values_Mode3     = [true false];
values_Mode4     = [true false];
%% Algorithm parameters
%%Variables
stepSize_R2 =    0.05;
stepSize_a  =     2.5;
values_R2   =   1.1: stepSize_R2  :1.5;
values_a31  =    10: stepSize_a   :20;
values_a32  =  47.5: stepSize_a   :57.5;         
values_a33  =  72.5: stepSize_a   :85;         
values_a41  = 112.5: stepSize_a   :132.5;
values_a42  = 127.5: stepSize_a   :140;
    

%%RF: Input for path calculation
paramRF.inputRange = [15 80];
paramRF.nTimeStepsRF = 16;      %number of timesteps for calculating RF paths

%%Evaluation paramaters
paramEval.threshold = tand(2); %tand(1);%0.05;
paramEval.minSectionLength = 4; %normal: 2
paramEval.nTimeStepsNew = 64;
  %check if nTimeStepsEval is suitable
  paramEval.factorTimeSteps = paramEval.nTimeStepsNew/paramRF.nTimeStepsRF;
  if round(paramEval.factorTimeSteps) ~= paramEval.factorTimeSteps
       error('The factor between nTimeStepsEval and nTimeStepsRF should be a whole number');
  end

%% STEP 1: Generate all valid sector angle combinations
% Find all geometrically valid paths for the proposed range/resolution
tic
validGeometries = getValidGeometries(values_a31,values_a32,...
                                values_a33,values_a41,values_a42);                      
toc
metaData.numberOfRFCalculations = size(validGeometries,2)*size(values_R2,2)*4;    %total number of calculations, for calculating progress
save(sprintf('Results/%s/%s_A_valid_geometries.mat',metaData.runID,metaData.runID),'validGeometries','metaData')

clearvars -except metaData validGeometries values_Mode3 values_Mode4 values_R2 paramRF paramEval R2Range
%% STEP 2&3: Calculate RF paths & Evaluate Feasibility
% Generate input
inputAngles = generateInputProfile(paramRF.inputRange,paramRF.nTimeStepsRF);

% counter
k = 1;
% Split range in parts, iterate

for Mode3 = values_Mode3
for Mode4 = values_Mode4
for R2 = values_R2(R2Range)
    %Generate save tag
    metaData.saveTag = sprintf('%s_M=%.0f%.0f_R2=%0.2f',metaData.runID,Mode3,Mode4,R2);
    %STEP 2: Generate RF paths
    [X_EE,variables] = generateRFPaths(R2,validGeometries,Mode3,Mode4,inputAngles);
    save(sprintf('Results/%s/%s_B_rigidly_foldable_paths.mat',metaData.runID,metaData.saveTag),'X_EE','variables','metaData')
    %STEP 3: Check Feasibility
    [feasible,flag] = evaluatePaths(X_EE,variables,paramEval,paramRF);
    save(sprintf('Results/%s/%s_C_feasiblePaths.mat',metaData.runID,metaData.saveTag),'feasible','metaData')
    %print progress % in command window
    sprintf('progress: %.1f %%',100*k/(size(values_R2(R2Range),2)*size(values_Mode3,2)*size(values_Mode4,2)))
    k = k+1;
end
end
end

%% Join results