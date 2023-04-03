%% Analysis: relation input speed output speed
clear
runID = '230201_08'; rbms =3;
load(sprintf('Results/%s/%s_D_analyzed.mat',runID,runID),'feasible','tbl','splittbl');

%% find min max horinode ranges

coordinates = feasible{2834}.X_EE;
%%
input = generateInputProfile([15 80],64);
speed_input1    = diff([input(1,:) input(1,1)]); 
speed_input2    = diff([input(2,:) input(2,1)]); 
t = 1:64;

for  i = 1:length(feasible)
    coordinates = feasible{i}.X_EE;
horistart(i) = feasible{i}.horiNodes(1);
horiend(i)   = feasible{i}.horiNodes(end);
speed_outputx        = diff([coordinates(:,1)' coordinates(1,1)]);
speed_outputz(i,:)   = diff([coordinates(:,2)' coordinates(1,2)]);
speed_output(i,:)    = vecnorm([speed_outputx;speed_outputz(i,:)]);
end

smallhori = [max(horistart) min(horiend)];
largehori = [min(horistart) max(horiend)];

avg_speed = mean(speed_output,1);
std_speed = std(speed_output);

%%
figure
set(gcf,'defaultAxesColorOrder',[0.35 0.4 .5;0.6 0 0]);

hold on

yyaxis left
p1 = plot(t,abs(speed_input1),'--','color',[0.7 0.8 1]);
p2 = plot(t,abs(speed_input2),'-.','color',[0.7 0.8 1]);
ylabel('Absolute input speed')
plot([smallhori(1),smallhori(1)],[0 3.5],'k-.');
plot([smallhori(2),smallhori(2)],[0 3.5],'k-.');
plot([largehori(1),largehori(1)],[0 3.5],'k--');
plot([largehori(2),largehori(2)],[0 3.5],'k--');
yyaxis right


% plot(t,speed_outputx);
% plot(t,speed_outputz);
% plot(t,speed_output);
e = errorbar(t,avg_speed,std_speed(1:2:end),'capsize',3);
ylabel('End effector speed (s_{EE})')

xlim([1,64])
title('Input angular speed vs. output speed')
% legend("d\rho_1/dt","\rho_2",'','','','',"s_{out}")
xlabel('t')
