%% filterIsomorphic()
% by Jim Sluijter
% Creates set of graphs to compare G1 with, calls function checkIsomorphic
% for every graph in this set

function isIso = filterIsomorphic(G1,G_SET)
    %Graphs are isomorphic if they exhibit an edge-preserving
    %vertex bijection
     %The set of vertex label Lv must coincide
     %i.e. the incomingNeighbours must be the same for all vertices
     
     if G1.g == 1           %G1 is the first of its generation, no comparison can be made 
        isIso = false;
        return
     else
        Set_nextgen = G_SET(G1.gen,:); %Create set of all comparable graphs in same generation
     end
     %Indices of these graphs
     V_g = find(~cellfun(@isempty,(Set_nextgen)));

     %default:
     isIso = false;
     %Check with ALL graphs in current writing generation
     for g = V_g
         G2 = Set_nextgen{g};
         tf =  checkIsomorphic(G1,G2);
         if tf
             isIso = true;
%              "graph ("+G1.gen+", "+G1.g+"): "+G1.label+"is isomorphic with graph ("+G2.gen+", "+G2.g+"): "+G2.label
             return
         end
     end
             