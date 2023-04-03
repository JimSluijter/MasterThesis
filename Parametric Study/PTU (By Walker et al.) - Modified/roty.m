%% License information

%This file is part of the matlabPTU framework.

% © 2021 ETH Zurich, Andreas Walker, D-MAVT, EDAC

%Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%%
function [matrix] = roty(ang)
%Returns a rotation matrix to rotate a vector by the specified angle (in degrees) around
%the y-axis
%Verifies that the angle is finite and calls rotateAboutArbitrary

% if not(isfinite(ang))
%     error('The specified angle is invalid. Expected a finite number.');
% end

% A=rotateAboutArbitrary(ang,[0;1;0]);
matrix = [cosd(ang) 0 sind(ang)
          0         1        0
         -sind(ang) 0 cosd(ang)];



end

