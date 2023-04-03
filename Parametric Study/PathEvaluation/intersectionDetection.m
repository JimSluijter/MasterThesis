%% intersection detection
% by Jim Sluijter
% Detects the majority of intersections by detecting if crease line c35
% crosses the fixed facet
function intersect = intersectionDetection(variables,input)
%extract nTimeSteps
nTimeSteps = size(input,2);
%Get variables
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
% allocate memory
X_EE = zeros(nTimeSteps,2);

% set default values
under = false;
above = false;
intersect = false;

for t = 1:nTimeSteps
    input13 = input(1,t);
    input23 = input(2,t);
    c35pos = getc35(Mode3,a31,a32,a33,a34,input13,input23);
    %Appoint values to memorize 
    if c35pos == -1
        under = true;
        if Mode3    %if Mode3 is true, this shouldn't happen unless there is an intersection
            intersect = true;
            return
        end
    elseif c35pos == 1
        above = true;
    end
    %If both have occured, there must be an intersection
    if under && above
        intersect = true;
        return
    end
end
        

%% Get c35;
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

end %function intersectionDetection
