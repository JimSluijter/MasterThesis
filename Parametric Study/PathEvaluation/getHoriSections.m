%% getHoriSections
%by Jim Sluijter
% This function extracts all horizontal sections. 
% A horizontal section is defined as a sequence of neighbouring horizontal 
% segments of at least minSectionLength long

% In the end, if multiple horizontal sections exist, they get split into
% multiple rows, so they can be analyzed separately


function horiSectionsI = getHoriSections(horiSegmentsV,minSectionLength)
nSegments = size(horiSegmentsV,1);
horiSectionsI = [];
%Test: filter out all sequences in horiSegmentsV that have fewer than
%minSectionLength in a row. 
%1. Add extra ones at the end if loop end in middle of sequence

%% Add extra segments to the end.
% In case a horizontal section starts halfway the loop, and crosses the
% start/end of the loop, we want the actual starting point to be
% considered, and not index 1. 
% Therefore, for further analysis:
    %if the vector ends with a horizontal segment, there is a possibility
    %for an section that crosses the loop
    %We check if this is the case, and if so, we add the extra entries to
    %the end of the vector:
        % e.g. a vector [1 1 0 0 0 1 1 1] becomes [1 1 0 0 0 1 1 1 (1) (1)]
        
if horiSegmentsV(end)  
    for i = 1:nSegments
        if horiSegmentsV(i)
            horiSegmentsV = [horiSegmentsV; true];
            horiSegmentsV(i) = false;               %makes sure the section is not repeated
        else
            horiSegmentsV = [horiSegmentsV; false]; %adds a false flag in orer to let the algorithm below know where to stop
            break %stop check
        end
    end
end


i = 1;
section = 0;

while i <= nSegments
    horiSectionV = zeros(size(horiSegmentsV)); %Allocate memory/reset values
    % We check for every node if it is the start of a sequence
    if horiSegmentsV(i)
        %if so: check if sequence is at least as long as minSectionLength
        longEnough = true;                  %initially set to true
        for j = 1:minSectionLength-1        %if any of the succeeding segments is not horizontal: sections is not long enough    
            if ~horiSegmentsV(i+j)          
                longEnough = false; 
                break
            end                   
        end
        %if so, write these sections onto vector horiSectionsV
        if longEnough
            j=0;
            section = section+1;
            while horiSegmentsV(i+j)
                horiSectionV(i+j) = true;
                j = j+1;
            end
            %Finding indices in right order:
            horiSectionsI{section} = find(horiSectionV);
            for ii = 1:length(horiSectionsI{section})
                if horiSectionsI{section}(ii) > nSegments
                    horiSectionsI{section}(ii) = horiSectionsI{section}(ii)-nSegments;
                end
            end
            
            i = i+j; %the rest of this particular section does not need to be checked again
        end
    end
    i = i+1;
end
end