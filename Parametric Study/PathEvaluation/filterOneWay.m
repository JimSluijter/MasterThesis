function oneWay = filterOneWay(Coordinates,horiNodesI)
% by Jim Sluijter
% checks if the horizontal section of the paths moves in a single direction
% if not, the path is discarded

horizontalX = Coordinates(horiNodesI,1);         %Extract X-values of nodes around horizontal sections with index g
    
% Check if one way:
%#1. sign(diff(horizontalX)) checks, for every value, the direction
 %if the sign is -1, the direction is in the negative direction and vice
 %versa
%#2. any(diff(...)) then checks if in this vector, there are any different
%values. If so, the horizontal section is not in a single direction
if any(diff(sign(diff(horizontalX))))
    oneWay=false;
else 
    oneWay = true;
end
end
