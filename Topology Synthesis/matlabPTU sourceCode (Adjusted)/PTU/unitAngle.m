%% License information

%This file is part of the matlabPTU framework.

% © 2021 ETH Zurich, Andreas Walker, D-MAVT, EDAC

%Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%%
function [unit,mat, gammas] = unitAngle(sector,dihedral)
%UNTITLED2 calculates unit angle of unit given by sector and dihedral
%angles
%   Input: sector: sector angles in counterclockwise order
%   dihedral: dihedral angles, counterclockwise order

if nargout>=2
    mat=eye(3);
end

if nargout==3
    gammas=zeros(1,2);
end

P1=[1;0;0];
P2=rotz(sector(end))*P1;
if nargout>=2
    mat=rotz(sector(end))*mat;
end
if nargout==3
    secondlast=P1;
end

for i=1:size(dihedral,2)
    if i==size(dihedral,2) && nargout==3
        second=P1;
        second=rotx(dihedral(end-i+1))*second; 
        second=rotz(sector(end-i))*second;
    end
    P2=rotx(dihedral(end-i+1))*P2; 
    P2=rotz(sector(end-i))*P2;
    if nargout==3
        secondlast=rotx(dihedral(end-i+1))*secondlast;
        secondlast=rotz(sector(end-i))*secondlast;
    end
    if nargout>=2
        mat=rotx(dihedral(end-i+1))*mat;
        mat=rotz(sector(end-i))*mat;
    end
end

if nargout==3 && size(dihedral,2)~=0
    gammas(1,1)=acosd(dot(second,P2));
    gammas(1,2)=acosd(dot(secondlast,P1));
end

unit=acosd(dot(P1,P2));

end

