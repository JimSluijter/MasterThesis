%% Categorize classes
% by Jim Sluijter
% This script is used to categorize the feasible paths for Rigid Body Mode
% FT, such that they are split in the three classes as indicated in the
% main paper. 

clear
runID = '230113_06'; %RBMs = 1:4
% runID = '230201_08';
load(sprintf('Results/%s/%s_D_analyzed.mat',runID,runID),'feasible','tbl','splittbl');
rbms = 3;
%% split path families
% Group I: 
    %horiSection is under z=0
    %There is a significant bend at the -x side of the horiSection:
        %look at acceleration in X-direction. That one is quite high
        %Also look at what moment it is. It is present at the moment the 
%Group II:
    %horisection is just below zero
    %smooth profile: acceleration is relatively constant

%Group III: 
    %horiSection is over z=0
    %We see that GroupIII is mostly defined by the sector angles in the first
    %index v3
ddx_threshold =0.02;
    
%allocate memory:
Group1 = boolean(zeros(size(feasible)));
Group2 = boolean(zeros(size(Group1)));
Group3 = boolean(zeros(size(Group1)));
for idx = splittbl{1,rbms}.index'
    
    horiNodesI = feasible{idx}.horiNodes;
    horiZsign= sign(feasible{idx}.X_EE(horiNodesI(1),2));
    
    %GROUP I
    [~,top] = max(feasible{idx}.X_EE(:,2));
    topX  = feasible{idx}.X_EE(top,1);
    horiX = [min(feasible{idx}.X_EE(horiNodesI,1)) max(feasible{idx}.X_EE(horiNodesI,1))];
    horiZ = min(feasible{idx}.X_EE(horiNodesI,2));
    if  horiZsign<0 && topX > horiX(2)
        Group1(idx) = true;
    end
    
    %GROUP II
    if  horiZsign<0 && topX <= horiX(2) && feasible{idx}.stepheight_body<0.5
        Group2(idx) = true;
    end
    
    %GROUP III
    if horiZsign>0
        Group3(idx) = true;
    end
    
    % check if a path is assigned to multiple groups:
    if Group1(idx)+Group2(idx)+Group3(idx) > 1
        error(" single path assigned to multiple groups; idx = "+num2str(idx) )
    end
    if Group1(idx)+Group2(idx)+Group3(idx) == 0
        warning("path idx "+num2str(idx)+"not assigned to group")
    end
end 

grp_tbl{1,rbms} = tbl(Group1,:);
grp_tbl{2,rbms} = tbl(Group2,:);
grp_tbl{3,rbms} = tbl(Group3,:);

%% heatmap
criteria_names = ["NumberOfFeasiblePaths"
                  "strideLength"
                  "stepHeight"
                  "contactRatio"
                  "error_Z"
                  "error_speed"
                  "Family" ];
              
makeHeatmaps(splittbl,rbms,"NumberOfFeasiblePaths");
for j = 1:3
    for k = 7;1:7
    makeHeatmaps(grp_tbl(j,:),rbms,criteria_names(k));
    sgtitle("Group "+num2str(j))
    end
end

%% compare dihedral angles
input = generateInputProfile([15 80],size(feasible{1}.X_EE,1));

%allocate memory
dihed34 = zeros(1,size(feasible{1}.X_EE,1));
dihed35 = zeros(size(dihed34));
dihed36 = zeros(size(dihed34));
dihed47 = zeros(size(dihed34));
dihed48 = zeros(size(dihed34));
dihed49 = zeros(size(dihed34));

dihed = zeros(size(feasible,2),size(input,2));
for i = [grp_tbl{1,rbms}.index' grp_tbl{2,rbms}.index' grp_tbl{3,rbms}.index'] %for all geometries in one of these families
    variables = feasible{i}.variables;
    R2  = variables(1);
    a31 = variables(2);
    a32 = variables(3);
    a33 = variables(4);
    a41 = variables(5);
    a42 = variables(6);
    Mode3 = variables(7);
    Mode4 = variables(8);
    a34 = 270 - (a31+a32+a33);
    a44 = 180 - a34;          
    a43 = 360 - (a41+a42+a44);

    for t = 1:length(input)
    [dihed36(t),dihed34(t),dihed35(t)] = PTU(Mode3,[a33],[],[a32],[],[a31, 90, a34],[input(1,t),input(2,t)]);
    [dihed47(t),dihed48(t),dihed49(t)] = PTU(Mode4,a42,[],[a41 a44],dihed34(t),a43,[]);
    end
    dihed(i,:) = dihed47; %only dihed47 now
end

dih47meanI   = mean(dihed(grp_tbl{1,rbms}.index',:));
stdI         = std(dihed(grp_tbl{1,rbms}.index',:));
dih47meanII  = mean(dihed(grp_tbl{2,rbms}.index',:));
stdII        = std(dihed(grp_tbl{2,rbms}.index',:));
dih47meanIII = mean(dihed(grp_tbl{3,rbms}.index',:));
stdIII       = std(dihed(grp_tbl{3,rbms}.index',:));
t = 1:64;

figure
hold on
plot(t,input(2,:))
errorbar([t t+64],[dih47meanI dih47meanI],[stdI stdI])
errorbar([t t+64],[dih47meanII dih47meanII],[stdII stdII])
errorbar([t t+64],[dih47meanIII dih47meanIII],[stdIII stdIII])