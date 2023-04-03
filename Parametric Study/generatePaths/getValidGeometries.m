%% getValidGeometries:
% by Jim Sluijter
% This function checks the geometric validity of every geometry by checking
% if it satisfies the conditions

function validGeometries = getValidGeometries(values_a31,values_a32,...
                                                values_a33,values_a41,values_a42)

% validGeometries = int8(zeros(5*10e7,6));  
validGeometries = zeros(5,10e7);
counter = 1;    

% k = 1;
for a31 = values_a31  
for a32 = values_a32
for a33 = values_a33
    a34 = 270 - (a31+a32+a33);
    if a34<=10 || a34>=160
       continue
    end
%     progress = k/12320
%     k = k+1;  
for a41 = values_a41
for a42 = values_a42
    %% Load variables & geometry checks
    %Check#1: Vertex 3 dimensions

    %Check#2: Vertex 4 dimensions
    a44 = 180 - a34;          
    a43 = 360 - (a41+a42+a44);
    if a43 <= 20 || a43 >= 160
        continue
    end
    
    % Save a variables
    validGeometries(:,counter) = [a31 a32 a33 a41 a42];
    counter = counter+1;
%     % Load last variable R2
%     nR2 = length(values_R2);
%     a_matrix = ones(size(values_R2'))*[a31 a32 a33 a41 a42];
%     validGeometries(counter:counter+nR2-1,:) = [values_R2' a_matrix];
%     counter = counter + nR2;
    
%     for R2 = values_R2
%         validGeometries(counter,:) = [R2 a31 a32 a33 a41 a42];
%         counter = counter+1;
%     end
end
end
end
end
end
validGeometries = validGeometries(:,1:counter-1);
