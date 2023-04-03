%% intersection detection
% by Jim Sluijter
% Filters out all paths where:
    % a dihedral angle is larger than 170 degrees (risk of
    % self-intersection)
    % all dihedral angles around a vertex are lower than 10 degrees (risk
    % of singularity)

function dihAngleStatus = filterDihedralAngles(variables,input)
%extract nTimeSteps
nTimeSteps = size(input,2);
%Get variables
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
% allocate memory
X_EE = zeros(nTimeSteps,2);

% set default values
dihAngleStatus = 0;

for t = 1:nTimeSteps
    input13 = input(1,t);
    input23 = input(2,t);
    %vertex3
    [dihed36,dihed34,dihed35] = PTU(Mode3,[a33],[],[a32],[],[a31, 90, a34],[input13 input23]);
    %We do not want any of the dihedral angles to be larger than 170;
    if any(abs([dihed34,dihed35,dihed36]) > 170)
        dihAngleStatus = 1;
        return
    end
    %vertex4
    [dihed47,dihed48,dihed49] = PTU(Mode4,a42,[],[a41 a44],dihed34,a43,[]);
    %We do not want that all dihedral angles around a vertex are small,
    %because that would be close to a singularity
    if all([dihed47,dihed48,dihed49] < 10)
        dihAngleStatus = -1;
        return
    %We also do not want any of the dihedral angles to be larger than 170
    %degrees
    elseif any(abs([dihed47,dihed48,dihed49]) > 170)
        dihAngleStatus = 1;
        return
    end
end