function pixelsperunit = getReferenceLength(testnr)
%this function is used to calculate the reference length in the picture by
%clicking two known parts of the mechanism.

%in this case, we take the length of R1 as a reference.
% R1 is 30 mm in the built model, and 1 in the computer model, so we measure the distance and divide by the 
image = imread(sprintf('Test%u-1.jpeg',testnr));
figure
imshow(image)
% Use ginput to capture the coordinates of the two points
disp('Select first point')
[x1, y1] = ginput(1);
disp('Select second point')
[x2, y2] = ginput(1);

vector = [x2-x1;y2-y1];
theta = acosd(dot(vector,[1;0])/(norm(vector)) );
disp(theta)
pixelsperunit = sqrt((x2-x1)^2 + (y2-y1)^2);