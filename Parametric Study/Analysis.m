%% Analysis
% by Jim Sluijter
% In this script, the remaining feasible paths are evaluated according to
% the defined criteria, in order to select the optimal geometry for the
% given design case

clear
runID = '230113_02';
runID = '230113_06';
runID = '230201_08';
% runID = '230201_09';
% runID = '230201_10';
% load(sprintf('Results/%s/%s_C_feasiblePaths_MERGED.mat',runID,runID)','feasible')
load(sprintf('Results/%s/%s_E_extra_requirements.mat',runID,runID)','feasible')
R1 = 1;
%% GRADE THE CRITERIA
for i = 1:length(feasible)
    if isstring(feasible{i})    %to filter out paths that were deemed unfeasible with higher resolution
        continue
    end
    coordinates = feasible{i}.X_EE;
    % In order to calculate the criteria correctly, we need to adjust the horiNodesI for the new resolution
    horiNodesI  = feasible{i}.horiNodes;
    threshold = 0.05;
    horiNodesI = extendHoriNodes(coordinates,horiNodesI,threshold);
    feasible{i}.horiNodes = horiNodesI;
    
    %% Properties
    % variables
    R2      = feasible{i}.variables(1);
    a31     = feasible{i}.variables(2);
    a32     = feasible{i}.variables(3);
    a33     = feasible{i}.variables(4);
    a41     = feasible{i}.variables(5);
    a42     = feasible{i}.variables(6);
    a34     = 270 - (a31+a32+a33);
    a44     = 180 - a34;          
    a43     = 360 - (a41+a42+a44);
    feasible{i}.a34 = a34;
    feasible{i}.a43 = a43;
    feasible{i}.a44 = a44;
    
    %size of mechanism
     %we take the max of either the size in X- or Y- direction, because
     %taking the size is mostly to outrule the extreme sizes.
    feasible{i}.XmechSize = 1+R1+feasible{i}.variables(1);
    feasible{i}.YmechSize = 1+abs(R1/tand(a34) + R2/tand(a43) );
    feasible{i}.mechSize = max(feasible{i}.XmechSize,feasible{i}.YmechSize);
    %% 1. Stride length
    % I can also define step size purely as the amount of horiSegments
    % compared to the total amount. 
    % The step size itself would be way less accurate, but the important
    % thing is: it shows the relative size between the 
    
    % For W_mech I could choose a couple things:
        %-just the X-size of the mechanism
        %-the Y-coordinate of v8, which determines the minimum size in
        %y-direction of the mechanism (when adheringto pure origami rules)
        
    firstHoriX = coordinates(horiNodesI(1),1);
    lastHoriX  = coordinates(horiNodesI(end),1);
    feasible{i}.stepSize = abs(firstHoriX-lastHoriX);

    feasible{i}.strideLength =  feasible{i}.stepSize/feasible{i}.mechSize;
    
    %% 2. Contact Ratio
    %-1 to only count segments instead of nodes
    feasible{i}.contactRatio = (length(feasible{i}.horiNodes) - 1)/length(feasible{i}.X_EE);
    %% 3. Error in speed
    Xcoords = coordinates(horiNodesI,1);
    dX = diff(Xcoords);

    feasible{i}.error_speed = abs( (max(dX)-min(dX)) / mean(dX));
    speedlabel = 'max. velocity difference';
    
    %% 4. Step height
    % we want a high point that is above the straight line section, so:
    Zcoords = feasible{i}.X_EE(:,2);
    %1. we find all nodes that have an X-coordinate between the firstHoriX
    %and lastHoriX
    bounds = sort([firstHoriX, lastHoriX]);
    lower = coordinates(:,1) < bounds(2);
    higher = coordinates(:,1) > bounds(1);
    inBetween = boolean(higher.*lower);
    inBetween(horiNodesI) = false;  %leave out the horiNodes
    if ~any(inBetween)
        feasible{i} = [];
        feasible{i} = "Invalid path";
        continue
    end
    %2. we find the maximum of their Z-values
    horiZValues = Zcoords(horiNodesI);
    feasible{i}.stepheight_path = max(max(abs(Zcoords(inBetween)-Zcoords(horiNodesI)')));
    feasible{i}.stepheight_body = min(abs(horiZValues));
    feasible{i}.stepHeight = min(feasible{i}.stepheight_path,feasible{i}.stepheight_body)/feasible{i}.mechSize;

    %% Error in Z
    % (normalized with stepSize)
    Zvals = coordinates(horiNodesI,2);
    feasible{i}.error_Z = abs(max(Zvals)-min(Zvals))/feasible{i}.stepSize;  
    
    %% Distance to 
end

%% filter out invalid paths
k=1;
for i = 1:length(feasible)
    if isstring(feasible{i})    %to filter out paths that were deemed unfeasible with higher resolution
        continue
    end
    feasible_new{k} = feasible{i};
    k = k+1;
end
clear feasible
feasible = feasible_new;
clear feasible_new

%% Extract to arrays
%criteria
stepHeight      = [];
strideLength    = [];
error_speed     = [];
error_Z         = [];
contactRatio    = [];
%sizes
YmechSize       = [];
XmechSize       = [];
mechSize        = [];
%step height
stepheight_path = [];
stepheight_body = [];
%variables
variables = [];
R2        = [];
a31       = [];
a32       = [];
a33       = [];
a41       = [];
a42       = [];
a34       = [];
a43       = [];
a44       = [];
Mode3     = [];
Mode4     = [];
rbm       = [];
index     = [];

%% Write in table
for i = 1:length(feasible)
    %check if still feasible
    if isstring(feasible{i})
        continue
    end

    %criteria
    strideLength    = cat(1,strideLength,   feasible{i}.strideLength);
    error_speed     = cat(1,error_speed,     feasible{i}.error_speed);
    stepHeight      = cat(1,stepHeight, feasible{i}.stepHeight);
    error_Z         = cat(1,error_Z,      feasible{i}.error_Z);
    contactRatio    = cat(1,contactRatio,    feasible{i}.contactRatio);
    %sizes
    YmechSize       = cat(1,YmechSize,            feasible{i}.YmechSize);
    XmechSize       = cat(1,XmechSize,            feasible{i}.XmechSize);
    mechSize        = cat(1, mechSize,            feasible{i}.mechSize);
    stepheight_path = cat(1,stepheight_path,      feasible{i}.stepheight_path);
    stepheight_body = cat(1,stepheight_body,      feasible{i}.stepheight_body);
    
    %variables
    index       = cat(1,index,i);
    variables   = cat(2,variables, feasible{i}.variables');
    R2          = cat(1, R2,       feasible{i}.variables(1));
    a31         = cat(1,a31,       feasible{i}.variables(2));
    a32         = cat(1,a32,       feasible{i}.variables(3));
    a33         = cat(1,a33,       feasible{i}.variables(4));
    a41         = cat(1,a41,       feasible{i}.variables(5));
    a42         = cat(1,a42,       feasible{i}.variables(6));
    a34         = cat(1,a34,       feasible{i}.a34);
    a43         = cat(1,a43,       feasible{i}.a43);
    a44         = cat(1,a44,       feasible{i}.a44);
    Mode3       = cat(1,Mode3,     feasible{i}.variables(7));
    Mode4       = cat(1,Mode4,     feasible{i}.variables(8));
    rbm         = cat(1,rbm,       string(feasible{i}.variables(7))+string(feasible{i}.variables(8)));
end

%Put in table
tbl   = table(index,R2,a31,a32,a33,a34,a41,a42,a43,a44,rbm,strideLength,error_speed,stepHeight,error_Z,...
    contactRatio,YmechSize,XmechSize,mechSize,stepheight_path,stepheight_body);
[groupidx,groupname] = findgroups(tbl.rbm);

% split for modes
NrOfRBMs = 4;
tags = ["11","10","01","00"];
for mode = 1:NrOfRBMs
splittbl{mode} = tbl(groupidx==find(strcmp(groupname,tags(mode))),:); %TT
end
% splittbl{2} = tbl(groupidx==find(strcmp(groupname,"10")),:); %TF
% splittbl{3} = tbl(groupidx==find(strcmp(groupname,"01")),:); %FT
% splittbl{4} = tbl(groupidx==find(strcmp(groupname,"00")),:); %FF

save(sprintf('Results/%s/%s_D_analyzed.mat',runID,runID),'feasible','tbl','splittbl');