%% License information
%This file is part of the matlabPTU framework.
% © 2021 ETH Zurich, Andreas Walker, D-MAVT, EDAC
%Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%%
function [dihedral1, dihedral2, dihedral3] = PTU(mode, sector1,...
    dihedralin1, sector2, dihedralin2, sector3, dihedralin3)
%Computes the unknown dihedral angles from the known sector and dihedral
%angles
%   Input: mode is a bool, true for "up", false for "down"
%   sector: array with sector angles, counterclockwise, in degrees
%   dihedralin: array with known dihedral angles, counterclockwise, in
%   degrees

%% Check #1: if sector is over 180 degrees
 %If correction is true for sector, the unit angle will be converted later
sums=[sum(sector1) sum(sector2) sum(sector3)];
correction=false(1,3);
for j=1:3
    correction(j)=(sums(j)>=180);
end
%% Check #2: Check sector dihedral
 %verifies that the dihedral vector contains one element
 %less than the sector angle
% checkSectorDihedral(sector1,dihedralin1);
% checkSectorDihedral(sector2,dihedralin2);
% checkSectorDihedral(sector3,dihedralin3);
    %% Calculate Unit Anlges
    unitAngles=zeros(3,1);
    [unitAngles(1), mat1, gammas1]=unitAngle(sector1,dihedralin1);
    [unitAngles(2), mat2, gammas2]=unitAngle(sector2,dihedralin2);
    [unitAngles(3), mat3, gammas3]=unitAngle(sector3,dihedralin3);
    %% Check #3: Check unit angles
    %Checks if (1) Sum of unit angles is not over 360 degrees
    %          (2) vertex is rigidly foldable (MIGHT CHANGE THIS to only occur
    %          at the max input angles if it takes too much time
    checkUnits(unitAngles);
%     if max(unitAngles)>median(unitAngles)+min(unitAngles)
%         error("Not rigidly foldable!");
%     end


%% Correct unit angles
for j=1:3
    if correction(j)
        unitAngles(j)=360-unitAngles(j);
    end
end

%% Calculate theta
theta1=sphericalCosines(unitAngles(1),unitAngles(2),unitAngles(3));
theta2=sphericalCosines(unitAngles(2),unitAngles(1),unitAngles(3));
theta3=sphericalCosines(unitAngles(3),unitAngles(1),unitAngles(2));

%% Assign signs and correct if needed
signs=zeros(3,2);
[signs(1,1), signs(1,2)]=signsFromB(mat1);
[signs(2,1), signs(2,2)]=signsFromB(mat2);
[signs(3,1), signs(3,2)]=signsFromB(mat3);
for j=1:3
    if correction(j)
        signs(j,:)=-signs(j,:);
    end
end

%% Calculate beta & delta
beta1=-signs(1,1)*sphericalCosines(gammas1(1,1),sector1(1,1),unitAngles(1));
beta2=-signs(2,1)*sphericalCosines(gammas2(1,1),sector2(1,1),unitAngles(2));
beta3=-signs(3,1)*sphericalCosines(gammas3(1,1),sector3(1,1),unitAngles(3));

delta1=-signs(1,2)*sphericalCosines(gammas1(1,2),sector1(1,end),unitAngles(1));
delta2=-signs(2,2)*sphericalCosines(gammas2(1,2),sector2(1,end),unitAngles(2));
delta3=-signs(3,2)*sphericalCosines(gammas3(1,2),sector3(1,end),unitAngles(3));

%% Calculate dihedral angles
if mode  
    dihedral1=beta3+180-theta1+delta2;
    dihedral2=beta1+180-theta2+delta3;
    dihedral3=beta2+180-theta3+delta1;
else
    dihedral1=beta3+theta1-180+delta2;
    dihedral2=beta1+theta2-180+delta3;
    dihedral3=beta2+theta3-180+delta1;
end
%% Check imaginary part of dihedral angle
%Checks if the imaginary part is non-negligible. If so, an error is
%produced, otherwise the imaginary part is set to 0.
dihedral1=checkImag(dihedral1);
dihedral2=checkImag(dihedral2);
dihedral3=checkImag(dihedral3);


%% Spherical cosines function 
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


end