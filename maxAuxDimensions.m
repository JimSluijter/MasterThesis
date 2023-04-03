%% Maximum auxiliary dimensions (FT)
% by Jim Sluijter
% This script is used to calculate the maximum auxiliary dimensions for
% Rigid Body Mode RBM_FT

%load
clear
load('C:\Users\jim\Documents\MATLAB\Master Thesis\Matlab code\V6 - Final Structure\Results\Final - Active Surface\feasible_2834.mat', 'feasible_2834')
variables = feasible_2834.variables;
Z_max = feasible_2834.stepheight_body;
Z4_max = sind(80)+Z_max;
input = generateInputProfile([15 80],64);
nTimeSteps = length(input);

%Variables
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

%calculate vectors
for t = 1:nTimeSteps
    rho13(t) = input(1,t);
    rho23(t) = input(2,t);

    %%PTU for both edges
    [rho36(t),rho34(t),rho35(t)] = PTU(Mode3,[a33],[],[a32],[],[a31, 90, a34],[rho13(t) rho23(t)]);
    [rho47(t),rho48(t),rho49(t)] = PTU(Mode4,[a42],[],[a41 a44],rho34(t),[a43],[]);

    %%Calculate edges c35 and c49
    %Known vectors
    c32 = [0;-1;0];
    N_z = [0;0;1];
    %calculate c_34
    c0_34 = rotateAboutArbitrary(a34,        N_z) *   c32;
    c34(:,t)   = rotateAboutArbitrary(rho23(t),       c32)  * c0_34;
    %calculate normal vector of facet A
    N_A(:,t)   = rotateAboutArbitrary(rho23(t),  [0;-1;0])  *   N_z;
    %calculate normal vector of facet D
    N_D(:,t) = rotateAboutArbitrary(rho34(t),     c34(:,t)) *   N_A(:,t);
    %calculate c35 and c49 by rotating over plane F_D
    c35(:,t)   = rotateAboutArbitrary(a33,        N_D(:,t)) *   c34(:,t);
    c49(:,t)   = rotateAboutArbitrary(-a41,       N_D(:,t)) *  -c34(:,t);
    
    %calculate c48
    N_B(:,t)   = rotateAboutArbitrary(rho47(t),   [0;-1;0]) * N_A(:,t);
    c48(:,t)   = rotateAboutArbitrary(a43,        N_B(:,t)) * [0;-1;0]; 
end

%find values of c35 and a49 at both steepest moments
% moment 1 is when c35 is at its steepest
[c35_zs,c35_idx] = min(c35(3,:));
c35steep = c35(:,c35_idx); 
L_c35max = Z_max/c35_zs;
% moment 2 is when c49 is at its steepest
[c49_zs,c49_idx] = min(c49(3,:));
c49steep = c49(:,c49_idx);
L_c49max = Z4_max/c49_zs;
% moment 3 is when c48 is at its steepest
[~,c48_idx] = min(c48(3,:));
c48steep = c48(:,c48_idx);

%for both creases, calculate the maximum angle for maximum shape (N_E_c35
%means the N_E at the moment c35 is at its steepest)

%around c35: E
N_E_c35 = rotateAboutArbitrary(rho35(c35_idx),  c35steep) * N_D(:,c35_idx);
intersection_E = cross(N_E_c35,N_z);
angle_E = acosd( dot(intersection_E,c35steep) / (norm(intersection_E)*norm(c35steep)) );
%aroujd c35: D
intersection_D35 = cross(N_D(:,c35_idx),N_z);
angle_D35 = acosd( dot(intersection_D35,c35steep) / (norm(intersection_D35)*norm(c35steep)) );

%around c49: D
intersection_D49 = cross(N_D(:,c49_idx),N_z);
angle_D49 = acosd( dot(intersection_D49,c49steep) / (norm(intersection_D49)*norm(c49steep)) );

%around c49: C
N_C_c49 = rotateAboutArbitrary(rho49(c49_idx) ,c49steep) * N_D(:,c49_idx);
intersection_C49 = cross(N_C_c49,N_z);
angle_C49 = acosd( dot(intersection_C49,c49steep) / (norm(intersection_C49)*norm(c49steep)) );

%around c48: C
N_C_c48 = rotateAboutArbitrary(rho48(c48_idx), c48steep) * N_B(:,c48_idx);
intersection_C48 = cross(N_C_c48,N_z);
angle_C48 = acosd( dot(intersection_C48,c48steep) / (norm(intersection_C48)*norm(c48steep)) );

clearvars -except R1 R2 a31 a32 a33 a34 a41 a42 a43 a44 L_c35max L_c49max angle_E angle_D35 angle_D49 angle_C49 angle_C48

%For the order of the cross product:
% if the plane is to the left of the crease: cross(N_z,N_p)
% if the plane is to the right of the crease: