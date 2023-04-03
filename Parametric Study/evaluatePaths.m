%% PathEvaluation
% Created by: Jim Sluijter, TU Delft
% This script loads the generated end-effector paths, and evaluates the
% following properties:
% #1. The script evaluates if a horizontal straight-line section is present
% in the path.
% #2. It evaluates if this straight-line no other nodes in the path exceed
% the Z-values of this straight-line section
% #3. It evaluates if this straight-line section is uni-directional
% #4. It evaluates if vertex 4 of the origami-crease pattern is does not
% exceed the Z-values of the straight-line section
function [feasible,flag] = evaluatePaths(X_EE,variables,paramEval,paramRF)
%% Set Parameters
threshold           = paramEval.threshold;
minSectionLength    = paramEval.minSectionLength;
nTimeSteps          = paramRF.nTimeStepsRF;
nTimeStepsNew       = paramEval.nTimeStepsNew;
factorTimeSteps     = paramEval.factorTimeSteps;

%% Recalculate input profile for intersection
inputAngles = generateInputProfile([15 80],nTimeSteps);

%% Allocate memory
% results = cell(size(X_EE,1),2);
feasible    = cell(size(X_EE));             %Used to seperately save feasible geometries
kk = 1;                                     %Index for 'feasible'
flag        = cell(size(X_EE));             %Used to save flags    
tic

%% Calculate the maximum of X_4
% This is calculated before the loop because it is independent of the
% geometry
Z_4max = sind(max(paramRF.inputRange));
%% Loop over every Rigidly Foldable Path
for k = 1:size(X_EE,1)
    %1. Calculate finite difference dZdX for every line segment
    coordinates = X_EE{k};    
    dZdX = diff([coordinates(:,2);coordinates(1,2)])./ diff([coordinates(:,1);coordinates(1,1)]);
        % The coordinate vectors are concatenated with their first entry to
        % make sure dZdX is calculated for the whole loop
        
    %2. Detect horizontal line sections
    horiSegmentsV = abs(dZdX)<threshold;    %indicates for every node if the succeeding segment is close to horizontal (true) or not (false)
    horiSegmentsI = find(horiSegmentsV);    %Indices of nodes that procede horizontal segments (i.e. low values of the finite differnce)
    %2a. check if any straight sections exist
        % this conditional is called before fiilterSectionLength to save
        % time
    if isempty(horiSegmentsI)
        flag{k,1} = "no horizontal segments present";
        continue
    end
    %2b. retrieve all sections with length => 'minSectionLength
    horiSectionsI = getHoriSections(horiSegmentsV,minSectionLength);
    
    
    for section = 1:numel(horiSectionsI)
        %3. return in case no sections with sufficient length exist
        if isempty(horiSectionsI{section})
            flag{k,section}  = "horizontal section(s) < minimum section length";
            continue    %go to next horizontal section
        end
        %4. Make sure horizontal section is higher than any other nodes in
        %path, and higher than the maximum Z of node 4.
        horiNodesI  = [horiSectionsI{section}; horiSectionsI{section}(end)+1];  %all nodes bordering horizontal segments in horiSections
            %in case this added node exceeds nSteps (this only happens when
            %it is the first node
            if horiNodesI(end) == nTimeSteps+1
                horiNodesI(end) = 1;
            elseif horiNodesI(end) > nTimeSteps+1
                error("Something went wrong")
            end
        horizontalZ = coordinates(horiNodesI,2)';       %Extract Z-values of nodes around horizontal sections with index
        %4a. make sure the horizontal section is either entirely above or
        %entirely below the mechanism
        if min(abs(horizontalZ)) < 0.05
            flag{k,section}  = "horiSection is too close to Z=0";
            continue
        end
        
        %4b. compare with vertex 4
        if min(horizontalZ)>0 && min(horizontalZ) < Z_4max
            flag{k,section}  = "Vertex 4 exceeds maximum Z-coordinates";
            continue    %go to next horizontal section
        end
        
        allOtherZ   = setdiff(coordinates(:,2),horizontalZ);  %A vector with all other Z-values to compare
        
        
        %4c. check if any of the following conditions occurs:
        if min(horizontalZ)>0 &&  min(horizontalZ) <= max(allOtherZ)
            flag{k,section}  = "straight section not on top";
            continue    %go to next horizontal section
        elseif max(horizontalZ)<0 && max(horizontalZ) >= min(allOtherZ)
            flag{k,section}  = "straight section not on bottom";
            continue    %go to next horizontal section
        elseif min(horizontalZ)<0 && max(horizontalZ)>0
            flag{k,section} = "straight section both above and below Z=0?";
            continue    %go to next horizontal section
        end
      
        %6. Make sure horizontal section is one-way
        oneWay = filterOneWay(coordinates,horiNodesI);
        if ~oneWay
            flag{k,section}  = "straight section not one way";
            continue    %go to next horizontal section
        end
        

        
        %7. SELF INTERSECTION FILTER
        % right now only look at the direction of c35
        intersect = intersectionDetection(variables{k},inputAngles);
        if intersect
            flag{k,section} = "self-intersection around c35";
            continue
        end
        
        %8. Outrule extreme dihedral angles
        %a. outrule singularities
%         b. outrule extreme dihedral angles
        dihAngleStatus = filterDihedralAngles(variables{k},inputAngles);
        if dihAngleStatus == 1
            flag{k,section} = "fold angle too large";
            continue
        elseif dihAngleStatus == -1
            flag{k,section} = "risk of singularity";
            continue          
        end
        
        % If we arrive at this point, the path has survived all checks. 
        % We save it in a new array in order to feed it to the next
        % function if needed, 
        flag{k,section}  = "Feasible";
        feasible{kk}.X_EE = X_EE{k};
        feasible{kk}.horiNodes = horiNodesI;
        feasible{kk}.variables = variables{k};
        feasible{kk}.sectionNr = section;
        kk = kk+1;
    end %section = 1:size(horiSectionsI,1)
end
feasible = feasible(1:kk-1,:); %adjust size of feasible
% info.time_evaluatePaths = toc;

% tic
%% Recalculate feasible mechanisms in higher resolution
if factorTimeSteps ~=1
    %we recalculate all paths with a highenumber of time steps
    
        %(nTimeStepsNew is defined at the top of the script)
    inputAngles = generateInputProfile(paramRF.inputRange,nTimeStepsNew);

    %Update info
    info.inputAngles = inputAngles;

    %Update X_EE and horiNodes
    for i = 1:length(feasible)
        %update X_EE
        try
            feasible{i}.X_EE = calcPath(feasible{i}.variables,inputAngles);
            %update horiNodes
            firstHori = feasible{i}.horiNodes(1);
            lastHori  = feasible{i}.horiNodes(end);
            firstHori = translateNodes(firstHori,factorTimeSteps);
            lastHori = translateNodes(lastHori,factorTimeSteps);
            %check if the horiSection wraps around the loop end
            if lastHori>=firstHori   % it doesn't
                feasible{i}.horiNodes = firstHori:lastHori';
            else % it does
                section1 = firstHori:nTimeStepsNew;
                section2 = 1:lastHori;
                feasible{i}.horiNodes = [section1 section2]';
            end  
        catch MExc
            feasible{i} = [];
            feasible{i} = "High resolution path not Rigidly foldable";
%             flag{k,section}  = "High resolution path not Rigidly foldable";
        end
    end%for i = 1:length(feasible)
end
% toc
% info.time_recalcPaths = toc;
% save(sprintf("Results/datasets/%s_feasiblePaths.mat",info.saveTag),'feasible','info');