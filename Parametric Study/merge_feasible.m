%% Merge feasible data
% This script is to merge all seperately stored data of one runID into a
% single data set to be used for evaluation of the paths.


runID = '230201_08_TT';
% values_R2 = 1.1:0.1:1.5;
values_R2 = 1.1:0.05:1.5;
% stepSize_R2 =   0.1;
% R2under1 = 0.2:stepSize_R2:1;
% values_R2   =  [R2under1, ones(1,length(R2under1))./fliplr(R2under1)];


temp_feasible = [];

for Mode3 = [true]
    for Mode4 = [true]
        for R2 = values_R2
            saveTag = sprintf('%s_M=%.0f%.0f_R2=%0.2f',runID,Mode3,Mode4,R2);

            load(sprintf('Results/%s/%s_C_feasiblePaths.mat',runID,saveTag))
            temp_feasible = cat(1,temp_feasible,feasible);
        end
    end
end

    
 feasible = temp_feasible;
 save(sprintf('Results/%s/%s_C_feasiblePaths_MERGED.mat',runID,runID)','feasible')
 
 save(sprintf('Results/%s/%s_C_feasiblePaths_MERGED.mat',runID,runID)','feasible')