%% Generate Input Profile
% By Jim Sluijter

% This code generates the input dihedral angles in creases c_13 and c_23
% This new version starts at a high input angle. This is because rigid
% foldability often fills with higher dihedral angles, and therefore, the
% algorithm spends less time generating paths where a part of the range is
% not rigidly foldable, creating a faster performing code.

function inputAngles = generateInputProfile(inputRange,nTimeSteps)

% Define phase shift
phaseShift = 90;
if round(nTimeSteps*phaseShift/360) ~= nTimeSteps*phaseShift/360
    warning(['Phase shift may be off due to innacurracy in rounding'
             'HINT: make sure nTimeStep*phaseShift/360 is an integer'])
end
phaseShift = round(nTimeSteps*phaseShift/360);

% Set input dihedral angle profiles
% for timesteps t from 1 to nTimeSteps
input13 = zeros(1,nTimeSteps);
for t = 1:nTimeSteps
    input23(t) = (inputRange(1)+inputRange(2))/2 + ...
        (inputRange(2)-inputRange(1))/2*cos( (t-1)*2*pi/nTimeSteps);
end
% Introduce phase shift in second input
% input 13 is 90 degrees ahead of input23

input13 = [input23(end-phaseShift+1:end)   input23(1:end-phaseShift)];

inputAngles = [input13;input23];
end

%% For trouble shooting: