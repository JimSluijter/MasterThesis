%% Image Processing code
% By Jim Sluijter
% This script is used to determine the end effector coordinates of the test
% results

clear
%% Inputs
%Test number
testnr =11;
%load PTU model coordinates
runID = '230201_08';
best_geometry_idx = 2834;
%% Extract coordinates from images
%Detect image scale
    % scale = getReferenceLength(testnr);
scale = 132.5;          %TEST11 (hardcoded to avoid repetition)
X_EE = zeros(16,2);     %Allocate memory
for i = 1:16
    %load image
    image = imread(sprintf('Test%u-%u.jpeg',testnr,i));                                                    
    %Detect green dots: find all pixels that have a minimum brightness for
    %green AND have an appropriate ratio of green vs. red and blue
    greenspots = image(:,:,2) > 50 &image(:,:,2) > 1.1*image(:,:,1) & image(:,:,2)> 1.1*image(:,:,3);
    
    %detect center coordinates
    [centers, radii] = imfindcircles(greenspots, [10 50], 'ObjectPolarity', 'bright', 'Sensitivity', 0.8);
    properties = regionprops(greenspots,'Area','Centroid');
    %filter out small green spots
    if length(properties) > 2 
        largespots = find(cat(2,properties(:).Area) > 200);
        properties = properties(largespots);
        if length(properties)>2
            figure
            imshow(imshow(greenspots))
            error('more than two spots detected')
        end
    end
    
    if i==5
        figure
        imshow(image)
        figure
        imshow(greenspots)
    end
    %detect right-most spot in order to find the right order of substraction
    centroids = cat(1,properties(:).Centroid);
    [~,maxidx] = max(centroids(:,1));
    minidx = mod(maxidx, 2) + 1;
    %Get coordinates of end-effector
    X_EE(i,:) = properties(maxidx).Centroid - properties(minidx).Centroid;
end

%Reverse the y-coordinate (image indexing of y-axis in matlab starts on top)
X_EE(:,2) = -X_EE(:,2);
%divide by image scale 
X_EE_experiment = X_EE./scale;

%% load PTU model coordinates for evaluation
load(sprintf('Results/%s/%s_D_analyzed.mat',runID,runID),'feasible','tbl','splittbl');
best_geometry = feasible{best_geometry_idx};
X_EE_model = best_geometry.X_EE;
X_EE_model_16 = X_EE_model(1:4:end,:);

%% Plot trajectories
figure(1)
clf
set(gcf,'defaultAxesColorOrder',[0.6 0 0; 0 0 1]);
subplot('Position',[0.1 0.4 0.8 0.7]) 
hold on
plot(X_EE_model(:,1),X_EE_model(:,2),'-','color',[0.6 0 0]);
plot(X_EE_model_16(:,1),X_EE_model_16(:,2),'o','color',[0.6 0 0]);
p3 = plot(X_EE_experiment(:,1),X_EE_experiment(:,2),'-ob');

%Appearance
title('Physical model vs. PTU model')
grid minor
grid on
axis equal
xlabel('X')
ylabel('Z')
axis([0.8 3 -1.2 0])
xticks([0.8,1.2,1.6,2.0,2.4,2.8])
yticks([-1.2,-0.8,-0.4,0])
%legend
p1 = plot(nan,nan,'-o','color',[0.6 0 0]);
legend([p1,p3],{'PTU model','Physical model'},'Location','northwest','Orientation','vertical');

% saveas(gcf, sprintf('%s_ExperimentvsPTU',runID), 'svg')
%% Plot distance between points
difference = X_EE_model_16 - X_EE_experiment;
for i = 1:16
error(i) = norm(difference(i,:))*sign(difference(i,1));
errorx(i) = difference(i,1);
errory(i) = difference(i,2);
end
p = polyfit(1:16,error,2);
x = 1:16;
fitline = p(1)*x.^2 + p(2).*x + p(3);
 
figure(1)
subplot('Position',[0.1 0.1 0.8 0.3]) 

hold on
yyaxis left
plot(1:16,X_EE_model_16(:,1)')
ylabel('XEE')
yyaxis right
% plot(1:16,error,'x');
% plot(1:16,fitline);
plot(1:16,errorx,'--*')
plot(1:16,errory,':*')
ylabel('error')
% plot(1:16,difference(:,1))
% plot(1:16,difference(:,2))
title('Correlation error vs. XEE')
legend('XEEmodel','errorx','errory','Location','northwest','Orientation','vertical');
xlabel('time t (s)')
% xticks([2,4,6,8,10,12,14,16])
xticklabels({'0','8','16','24','32','40','48','56','64'})
