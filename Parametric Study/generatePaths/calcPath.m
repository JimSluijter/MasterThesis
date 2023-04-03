%% function to calculate entire path
function X_EE = calcPath(variables,input)
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

X_EE = zeros(nTimeSteps,2);
for t = 1:nTimeSteps
    input13 = input(1,t);
    input23 = input(2,t);
    [~,X_EE(t,:)] = CalculateX_EE(Mode3,Mode4,R2,a31,a32,a33,a34,a41,a42,a43,a44,input13,input23);
end
