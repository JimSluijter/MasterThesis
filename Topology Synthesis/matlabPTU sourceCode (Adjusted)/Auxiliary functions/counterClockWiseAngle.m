%% License information

%This file is part of the matlabPTU framework.

% © 2021 ETH Zurich, Andreas Walker, D-MAVT, EDAC

%Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%%
function [ang] = counterClockWiseAngle(vec1,vec2,normalVec)
%COUNTERCLOCKWISEANGLE returns angle between vec1 and vec2, measured in
%counterclockwise direction. Accordingly, the possible output ranges from 0
%to 360. The last input argument specifies the normal plane on which the
%two vectors lie (and therefore the direction in which the angle is
%measured)
%   Output in degrees.
%   INPUT
%   vec1: 3-by-1 vector
%   vec2: 3-by 1 vector
%   Problem Sketch:
%
%   0------------> vec2
%   |phi/
%   |--/
%   |
%   |
%   \/
%   vec1

if (norm(vec1)<1e-9) || (norm(vec2)<1e-9)
    error('At least one vector is too short for angle calculation!');
end

if norm(vec1-vec2,2)<1e-12
    %the two vectors are identical. The counterclockwise angle inbetween is
    %either 0 or 360 (equivalent). We return 360.
    ang=360;
    return;
end

auxVec=rotateAboutArbitrary(90,normalVec)*vec1;

%capping the argument to not go below -1 or above 1
normalAngle=acosd(max([-1,min([dot(vec1,vec2)/(norm(vec1)*norm(vec2)),1])]));

if dot(auxVec,vec2)>0
    ang=normalAngle;
else
    ang=360-normalAngle;
end

end

