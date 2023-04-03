function horiNodesI = extendHoriNodes(coordinates,horiNodesI,threshold)
%by Jim Sluijter
% In this function, the nodes next to horiNodesI are examined to see if
% they are horizontal as well

%nSteps
nSteps = length(coordinates);

%% extend on the positive side
extending = true;
%start with last horiNode
    idx = horiNodesI(end);
while extending
    
    %if idx wraps around
    if idx > nSteps
        idx = idx-nSteps;
    end
    next_idx = idx+1;
    if next_idx > nSteps
        next_idx = next_idx-nSteps;
    end
        
    %calculate finite difference with the next node (dZdX is a scalar)
    dZdX = diff([coordinates(idx,2);coordinates(next_idx,2)]) ./ diff([coordinates(idx,1);coordinates(next_idx,1)]);
    if abs(dZdX)<threshold
        %horiNodesI should be extended!
        horiNodesI = [horiNodesI, next_idx];
        idx = idx+1; %move on to next idx
    else
        %horiNodesI should not be extended
        extending = false;
    end
end

%% extend on negative side
extending = true;
%start with first horiNode
idx = horiNodesI(1);
while extending
    
    %if idx wraps around
    if idx < 1
        idx = idx+nSteps;
    end
    next_idx = idx-1;
    if next_idx < 1
        next_idx = next_idx+nSteps;
    end
    
    %calculate finite difference with the next node (dZdX is a scalar)
    %(note that next_idx and idx are switched here)
    dZdX = diff([coordinates(next_idx,2);coordinates(idx,2)]) ./ diff([coordinates(next_idx,1);coordinates(idx,1)]);
    if abs(dZdX)<threshold
        %horiNodesI should be extended!
        horiNodesI = [next_idx horiNodesI];
        idx = idx-1; %move on to next i
    else
        %horiNodesI should not be extended
        extending = false;
    end
end