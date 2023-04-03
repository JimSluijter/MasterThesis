%% License information

%This file is part of the matlabPTU framework.

% © 2021 ETH Zurich, Andreas Walker, D-MAVT, EDAC

%Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%% 
function [theta] = sphericalCosines(opposite, side1, side2)
%SPHERICALCOSINES Uses Law of spherical cosines to calculate angle between
%two units from angles in unit triangle
theta=real(acosd((cosd(opposite)-cosd(side1)*cosd(side2))/(sind(side1)*sind(side2))));

if abs(imag(acosd((cosd(opposite)-cosd(side1)*cosd(side2))/(sind(side1)*sind(side2)))))>1e-2
    fprintf('something gone wrong. opposite=%e side1=%e side2=%e argument=%e\n',opposite,side1,side2,...
        (cosd(opposite)-cosd(side1)*cosd(side2))/(sind(side1)*sind(side2)));
    error('invalid spherical triangle');
end

if opposite>=180
    theta=360-theta;
end
end

