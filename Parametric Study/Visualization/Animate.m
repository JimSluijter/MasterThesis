%% Visualise mechanism
%by Jim Sluijter
% This script uses the adjusted method animateEEpath in the origami.m class
% to animate the chosen geometry
clear o

%three groups:
% runID = '230201_08'; i =  2834; %GROUPI
% runID = '230113_06'; i = 22586; %GROUPII
runID = '230113_06'; i = 22427; %GROUPIII

% update origami file line 2593 to compensate for bug!

  %% load parameters
load(sprintf('Results/%s/%s_D_analyzed.mat',runID,runID),'feasible','tbl','splittbl')
input = generateInputProfile([15 80],16);
%
R2      = feasible{i}.variables(1);
a31     = feasible{i} .variables(2);
a32     = feasible{i}.variables(3);
a33     = feasible{i}.variables(4);
a41     = feasible{i}.variables(5);
a42     = feasible{i}.variables(6);
Mode3   = feasible{i}.variables(7);
Mode4   = feasible{i}.variables(8);

%% Load crease pattern object
load('Results/Parallel_Stencil.mat')
%Update crease pattern object
o.animationView = [170 15];

o.vertexToGripWith = 8;
o.setModeTo(3,Mode3); o.setModeTo(4,Mode4);
o.updateDimensions('R2',R2,'a31',a31,'a32',a32,'a33',a33,'a41',a41,'a42',a42);%true false
%Update nodes
o.vertices(2).position0(2) = o.vertices(7).position0(2);
o.vertices(18).position0(2) = o.vertices(7).position0(2);
o.resetAdjacentAuxiliaries(7);
o.vertices(22).position0(2) = o.vertices(7).position0(2);
o.vertices(22).position0(1) = o.vertices(8).position0(1);
%22427
if i == 22427
    L36 = 0.35; L35 = 0.5; L49 = 0.5;
    a34 = 270 - (a31+a32+a33);
    a44 = 180 - a34;          
    a43 = 360 - (a41+a42+a44);
    o.vertices(6).position0 = [-L36*cosd(a31);L36*sind(a31);0];
    o.vertices(5).position0 = [-L35*cosd(a31+a32);L35*sind(a31+a32);0];
    o.vertices(9).position0 = [o.vertices(4).position0(1)+L49*sind(a43+a42)
                               o.vertices(4).position0(2)-L49*cosd(a43+a42);0];
    o.resetAdjacentAuxiliaries(5);
    o.resetAdjacentAuxiliaries(6);
    o.resetAdjacentAuxiliaries(9);
end

%% Plots: still 
% o.plotOrigamiGraph();
% o.plotOrigamiUndeformed();title(num2str(i));

o.calculateSectorLengths();
% o.setDihedralAngles([1 1]);
% o.determineOrigami();
% o.plotOrigamiDeformed()
%% plots: animation
o.animateEEpath([15 80],[15 80],sprintf('Results/figures/animation_%s/animation_%s_i=%u_B.gif',runID,runID,i),[-2 3.5],[-3.5 3],[-1.5,1.5])
title(num2str(i))