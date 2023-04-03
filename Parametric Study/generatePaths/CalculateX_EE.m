%% Calculate EE cooordinates for single configuration
% 24/11/22
% See notes for 24-11 for further description
% This is a further simplified version that calculates the EE coordinates,
% but omits the y-coordinate since this is not important for the parametric
% sweep. Therefore a simple 2D rotation matrix can be used
% X_4 is also calculated, since it is needed in the filter
% Try-catch scheme is deleted. It now just gives an error, which is picked
% up by the try-catch in ParametricSweep.m

function [X_4,X_EE] = CalculateX_EE(Mode3,Mode4,R2,a31,a32,a33,a34,a41,a42,a43,a44,input13,input23)
%% PTU around vertex 3 and 4
% PTU at v3 calculates dihed34, which is an input at PTU at v4
[dihed36,dihed34,dihed35] = PTU(Mode3,a33,[],a32,[],[a31, 90, a34],[input13 input23]);
    %We do not want any of the dihedral angles to be larger than 170;
    if any(abs([dihed34,dihed35,dihed36]) > 175)
        error('fold angle too large in vertex 3')
    end
%Vertex 4
[dihed47,dihed48,dihed49] = PTU(Mode4,a42,[],[a41 a44],dihed34,a43,[]);
    if all([dihed47,dihed48,dihed49] < 5)
        error('risk of singularity in vertex 4')
    elseif any(abs([dihed47,dihed48,dihed49]) > 175)
        error('fold angle too large in vertex 4')
    end
%% Coordinates vertex 4 and EE
X_4 = [cosd(input23);sind(input23)]; %R1*[cos(input23);sin(input13)] but R1=1 always
X_EE = X_4 + R2*[cosd(input23+dihed47);sind(input23+dihed47)];
