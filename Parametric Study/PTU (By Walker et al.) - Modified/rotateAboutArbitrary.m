%% License information

%This file is part of the matlabPTU framework.

% © 2021 ETH Zurich, Andreas Walker, D-MAVT, EDAC

%Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%%
function [R] = rotateAboutArbitrary(theta, ax)
%rotateAboutArbitrary returns matrix to perform rotations around angle
%theta around arbitrary axis ax.
ax=ax/norm(ax);
u=ax(1,1);
v=ax(2,1);
w=ax(3,1);
R=[cosd(theta)+u^2*(1-cosd(theta)) u*v*(1-cosd(theta))-w*sind(theta) u*w*(1-cosd(theta))+v*sind(theta);...
    u*v*(1-cosd(theta))+w*sind(theta) cosd(theta)+v^2*(1-cosd(theta)) v*w*(1-cosd(theta))-u*sind(theta);...
    u*w*(1-cosd(theta))-v*sind(theta) v*w*(1-cosd(theta))+u*sind(theta) cosd(theta)+w^2*(1-cosd(theta))];
end

