%% Check how many straight sections exist
% by Jim Sluijter
% This function splits the geometry into two different geometry instances
% in case there are multiple horizontal sections, such that each horizontal
% section can be analyzed seperately
function horiSectionsNew = splitMultipleSections(horiSectionsI)

if isempty(horiSectionsI)
    horiSectionsNew = [];
    return
end

difference = diff([horiSectionsI;horiSectionsI(1)]); %get difference
%Takes the number of times a section ends (= number of sections)

sectionEndI = find(difference~=1); %indicates the indices of horiSections where the sections end
nSections = sum(difference~=1); %simple as that:) one line

% Check if nSections > 0:
%if nSections == 0, we can state that there is not horizontal section   
%split sections in two solutions :)
horiSectionsNew = cell(nSections,1);
%for the first
horiSectionsNew{1} = horiSectionsI(1:sectionEndI(1));
%for all following
for i = 2:nSections
    horiSectionsNew{i} = horiSectionsI(sectionEndI(i-1)+1:sectionEndI(i));
end
end



        
   