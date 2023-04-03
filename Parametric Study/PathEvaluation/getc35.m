%% Get c35;
%by Jim Sluijter
% function to find the position of the end vertex of crease c35 when crease
% length is of unit length
function c35pos = getc35(Mode3,a31,a32,a33,a34,rho13,rho23)
%Use ptu to get rho34
[~,rho34,~] = PTU(Mode3,[a33],[],[a32],[],[a31, 90, a34],[rho13 rho23]);
%Known vectors
c32 = [0;-1;0];
N_z = [0;0;1];
%calculate c_34
c0_34 = rotateAboutArbitrary(a34,        N_z) *   c32;
c34   = rotateAboutArbitrary(rho23,      c32) * c0_34;
%calculate normal vector of facet A
N_A   = rotateAboutArbitrary(rho23, [0;-1;0]) *   N_z;
%calculate normal vector of facet E
N_E   = rotateAboutArbitrary(rho34,      c34) *   N_A;
%calculate c35
c35   = rotateAboutArbitrary(a33,        N_E) *   c34;

%Check position of c35 is neither, above, or below facet A
if c35(1) >= 0 || c35(2) >= 0
    c35pos =  0;    
elseif c35(1) < 0 && c35(2) < 0 && c35(3) > 0
    c35pos =  1;
elseif c35(1) < 0 && c35(2) < 0 && c35(3) <= 0
    c35pos = -1;
else
    
end %if statement

end %function getc35