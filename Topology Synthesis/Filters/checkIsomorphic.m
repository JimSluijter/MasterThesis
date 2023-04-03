%% checkIsomorphic()
% by Jim Sluijter

% Checks for two graphs, G1 and G2, if they are Isomorphic   

function isIso = checkIsomorphic(G1,G2)

%% early return if no match is possible
%checks number of nodes, number of edges
NnodesG1 = G1.numberOfVertices('omitDeleted');
NnodesG2 = G2.numberOfVertices('omitDeleted');

NedgesG1 = G1.numberOfEdges('omitDeleted');
NedgesG2 = G2.numberOfEdges('omitDeleted');

if NnodesG1 ~=NnodesG2 || NedgesG1 ~= NedgesG2
    isIso = false;
    return;
end

%% 
%extract node indices of non-deleted nodes
NodesG1 = find(~cat(2,G1.vertices.deleted));

NodesG2 = find(~cat(2,G2.vertices.deleted));
%% create unique index for every combination of predecessors
%first we convert all combs of predecessors to strings
%in case there are no predecessors ==> '0'
predG1 = [];
predG2 = [];
for a = NodesG1 
    p = G1.vertices(a).incomingNeighbours;
    if isempty(p)
        pstr = "0";
    else 
        pstr = num2str(p);
    end
    predG1 = [predG1;pstr];
end
for a = NodesG2 
    p = G2.vertices(a).incomingNeighbours;
    if isempty(p)
        pstr = "0";
    else 
        pstr = num2str(p);
    end
    predG2 = [predG2;pstr];
end
%now we can create unique indices for every predecessor combo
propTable = [predG1;predG2];
[~,~,nodeCat] = unique(propTable);
%a~d put them in two vectors next to eachother
nodeCat = uint64(reshape(nodeCat, [size(predG1,1) 2]));

% Check if node index distributions are equal; otherNowise, there is no
 % possible permutation that matches the node properties
    matchPossible = isequal(accumarray(nodeCat(:, 1), 1), accumarray(nodeCat(:, 2), 1));

%  if ~matchPossible 
%      isIso = false;
%      return
%  end

isIso = matchPossible; 