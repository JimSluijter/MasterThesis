%% checkSymmetric()
% Function to check symmetry between graphs G1 and G2
% Created by Jim Sluijter

function isSym = checkSymmetric(G1,G2)

%% early return if no match possible
%checks number of nodes, number of edges
NnodesG1 = G1.numberOfVertices('omitDeleted');
NnodesG2 = G2.numberOfVertices('omitDeleted');

NedgesG1 = G1.numberOfEdges('omitDeleted');
NedgesG2 = G2.numberOfEdges('omitDeleted');

if NnodesG1 ~=NnodesG2 || NedgesG1 ~= NedgesG2
    isSym = false;
    return;
end

%% extract node positions
NodesG1 = find(~cat(2,G1.vertices.deleted)); %selects only the non-deleted nodes
NodesG2 = find(~cat(2,G2.vertices.deleted));
G1pos = cat(2,G1.vertices(NodesG1).position0);
G1pos = G1pos(1:2,:);
G2pos = cat(2,G2.vertices(NodesG2).position0);
G2pos = G2pos(1:2,:);

%for all vertices in G2, a column vector must be present in G1 that has the
%same coordinates, but the sign of y is switched
G1pos_mirror = G1pos;
G1pos_mirror(2,:) = -G1pos(2,:);


MAP = zeros(1,size(G2pos,2));
for a = 1:size(G2pos,2)
    map = find(all(G2pos(:,a) == G1pos_mirror));
    if isempty(map)
        isSym = false;
        return
    end
    MAP(a) = map;
end

if all(sort(MAP) == 1:size(G2pos,2))
    isSym = true;
else
    error("check checkSymmetric")
end