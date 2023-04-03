%% Scatter plot
clear
% runID = '230113_05'; %RBMs = 3;
runID = '230201_08'; rbms =3; 1:4;
% runID = '230201_10'; rbms =3;
load(sprintf('Results/%s/%s_D_analyzed.mat',runID,runID),'feasible','tbl','splittbl');


%%
errorMargin = 0.01; %(is one percent)

str_mode = ["TT","TF","FT","FF"];
str_vari = ["R2","a31","a32","a33","a41","a42"];
str_vari_text = ["R_{2}","\alpha_{31}","\alpha_{32}","\alpha_{33}","\alpha_{41}","\alpha_{42}"];

%labels
label_error = '\epsilon_z';'normalized error \epsilon [-]';
label_height= 'normalized step height h_{step}';
label_size  = 'normalized step size s_{step}';
label_speed = '\epsilon_v';'maximum speed difference in step section';

i_best = 4641; %2834
%% STEP SIZE VS. STEP HEIGHT (FOR WALKING)
f0 = figure;
clf
% sgtitle("\epsilon_Z vs. \epsilon_s")
f0.Color    = 'w';
for mode = rbms
%Extract feasible paths with small error 
nosmallsector = splittbl{mode};%(splittbl{mode}.a31 > 12.5,:);
selected = splittbl{mode}(splittbl{mode}.contactRatio > 0.25 & splittbl{mode}.a31 > 12.5,:);%splittbl{mode}(splittbl{mode}.contactRatio > 0.25,:);%splittbl{mode}(splittbl{mode}.error_Z < errorMargin,:);
% large_propPort = splittbl{mode}(splittbl{mode}.propulsionPortion > 0.3,:);

%get direction of v48:
av48 = selected.a43- 90;

%Scatter plot
% subplot(2,2,mode)
title('L_s');
hold on
grid on
% s1{mode} = scatter(splittbl{mode}.strideLength,   splittbl{mode}.stepHeight,[],                      'filled');%splittbl{mode}.propulsionPortion,'filled');%log(splittbl{mode}.normalizedError),'filled');%splittbl{mode}.a42,'filled')
% s2{mode} = scatter(selected.strideLength,               selected.stepHeight,25,'filled');%large_propPort.propulsionPortion,'filled');%log(large_propPort.normalizedError),'filled');%selected.a42,'filled')%

s1{mode} = scatter(nosmallsector.error_Z, nosmallsector.error_speed,[],                 'filled');%splittbl{mode}.propulsionPortion,'filled');%log(splittbl{mode}.normalizedError),'filled');%splittbl{mode}.a42,'filled')
s2{mode} = scatter(selected.error_Z,      selected.error_speed,25,selected.strideLength,'filled');%large_propPort.propulsionPortion,'filled');%log(large_propPort.normalizedError),'filled');%selected.a42,'filled')%
colormap(turbo);
c=colorbar;
%     c.Title = label_speed;
%     c.Limits = [min(tbl.error_speed) max(tbl.error_speed)];
    %Appearance of INFEASIBLE data points
s1{mode}.MarkerEdgeColor = 'k';
%     s1{mode}.MarkerFaceColor =  clr(mode);
s1{mode}.MarkerEdgeAlpha = 0.1;
s1{mode}.MarkerFaceAlpha = 0.005;
s1{mode}.Marker = '+';
s1dtRow = [dataTipTextRow("\DeltaX_{step}/w_{mech}",splittbl{mode}.error_Z)
           dataTipTextRow("h_{step}/w_{mech}",      splittbl{mode}.error_speed)
           dataTipTextRow("index",                  splittbl{mode}.index)
           dataTipTextRow("R2",                     splittbl{mode}.R2)
           dataTipTextRow("a31",                    splittbl{mode}.a31)
           dataTipTextRow("a32",                    splittbl{mode}.a32)
           dataTipTextRow("a33",                    splittbl{mode}.a33)
           dataTipTextRow("a41",                    splittbl{mode}.a41)
           dataTipTextRow("a42",                    splittbl{mode}.a42)
           dataTipTextRow("RBM",                    splittbl{mode}.rbm)];
s1{mode}.DataTipTemplate.DataTipRows = s1dtRow;
%Appearance of FEASIBLE data points
s2{mode}.MarkerEdgeColor = 'k';
% s2{mode}.MarkerFaceColor = colormap;
s2{mode}.MarkerEdgeAlpha = 0.8;
s2{mode}.MarkerFaceAlpha = 0.8;
s2{mode}.LineWidth = 0.5;
s2dtRow = [dataTipTextRow("e_Z",                    selected.error_Z)
           dataTipTextRow("e_s",                    selected.error_speed)
           dataTipTextRow("L_s",                    selected.strideLength)
           dataTipTextRow("C",                      selected.contactRatio)
           dataTipTextRow("index",                  selected.index)
           dataTipTextRow("R2",                     selected.R2)
           dataTipTextRow("a31",                    selected.a31)
           dataTipTextRow("a32",                    selected.a32)
           dataTipTextRow("a33",                    selected.a33)
           dataTipTextRow("a41",                    selected.a41)
           dataTipTextRow("a42",                    selected.a42)];  
s2{mode}.DataTipTemplate.DataTipRows = s2dtRow;   
% axis('equal')
% xlim([0 0.5]);
% ylim([0 0.3]);
xlim([0,0.03]);
ylim([0.3,1.2]);
xlabel(label_error)
ylabel(label_speed)
 
s3{mode} = scatter(tbl.error_Z(i_best),tbl.error_speed(i_best),100);
    %Best
    
    %get color of best data point
    scatterChildren = get(gca, 'Children');
    bestidx = find(selected.index == i_best);
    maxcolorval = max(scatterChildren(2).CData);
    mincolorval = min(scatterChildren(2).CData);
    colorval = scatterChildren(2).CData(bestidx);
    coloridx = floor( 256*(colorval-mincolorval)/(maxcolorval-mincolorval));
    colormap = turbo(256); 
    s3{mode}.MarkerFaceColor = colormap(coloridx,:);
    s3{mode}.Marker = 'hexagram';
    s3{mode}.MarkerEdgeColor = 'k';
    s3{mode}.LineWidth = 1;
    s3{mode}.MarkerEdgeAlpha = 1.0;
    s3{mode}.MarkerFaceAlpha = 1.0;

%     
%     %Best
%     s3{mode}.Marker = 'hexagram';
%     s3{mode}.MarkerEdgeColor = 'k';
%     s3{mode}.MarkerFaceColor = 'c';
%     s3{mode}.MarkerEdgeAlpha = 1.0;
%     s3{mode}.MarkerFaceAlpha = 1.0;
%         s3dtRow = [dataTipTextRow("\DeltaX_{step}/w_{mech}",selected.strideLength)
%                dataTipTextRow("h_{step}/w_{mech}",      selected.error_speed)
%                dataTipTextRow("index",                  selected.index)
%                dataTipTextRow("R2",                     selected.R2)
%                dataTipTextRow("a31",                    selected.a31)
%                dataTipTextRow("a32",                    selected.a32)
%                dataTipTextRow("a33",                    selected.a33)
%                dataTipTextRow("a41",                    selected.a41)
%                dataTipTextRow("a42",                    selected.a42)
%                dataTipTextRow("RBM",                    selected.rbm)];  
%     s3{mode}.DataTipTemplate.DataTipRows = s2dtRow;   
%Labels

end
legend(["C \leq 0.25"; "C > 0.25";"best option"])
text(0.001,0.48,"G_{4641}")
%%
figure
clf
hold on

for i = 22427;[2834];[4105 6040]
plot([feasible{i}.X_EE(:,1);feasible{i}.X_EE(1,1)],[feasible{i}.X_EE(:,2);feasible{i}.X_EE(1,2)],'o-')
title("index:"+num2str(i))
axis('equal')
subtitle(sprintf("R_{2}=%.1f deg | a_{31}=%.0f | a_{32}=%.1f | a_{33}=%.1f | a_{41}=%.1f | a_{42}=%.1f",feasible{i}.variables(1:6))+sprintf("RBM:%u%u",feasible{i}.variables(7:8)))

end
