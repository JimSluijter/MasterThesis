%% Function to make heatmaps

function [f] = makeHeatmaps(splittbl,rbms,criterium,runID)
% Inputs:
    %splittbl:  structure of tables, one for each mode
    %rbms:      which RBMs to plot a heatmap for
    %criterium: string denoting the criterium for which the heatmap is
    %plotted
    %runID:     (optional), to save figures
    
%% Setup   


    
% Strings
str_mode = ["TT","TF","FT","FF"];
str_vari = ["R2","a31","a32","a33","a41","a42"];
str_vari_text = ["R_{2}","\alpha_{31}","\alpha_{32}","\alpha_{33}","\alpha_{41}","\alpha_{42}"];

%decide criteria index
criteria_names = ["NumberOfFeasiblePaths"
                  "strideLength"
                  "stepHeight"
                  "contactRatio"
                  "error_Z"
                  "error_speed"
                  "Family" ];
criteria_idx = find(criteria_names == criterium);
%call error if input criterium is something else
if isempty(criteria_idx)
    error('criterium not known')
end

figure_titles = ["Number of feasible paths"
                 "Max. stride length"
                 "Max. step height"
                 "Max. contact ratio" 
                 "Min. height error"
                 "Min. speed error"
                 "Family"];
            
heatmap_titles = ["N_{paths}";"l_s";"h_s";"C_{ }";"log(\epsilon_{z})";"log(\epsilon_{s})";"Family"];


colormethod = ["nothing","max","max","max","min","min","nothing"];

%Show all data             
XData = string(0.6:0.1:3);
YData = string(10:5:170);
%Xticks and Yticks
XLabels = repmat("",size(XData));
XLabels(1:2:end) = XData(1:2:end);
YLabels = repmat("",size(YData));
YLabels(1:4:end) = YData(1:4:end);

colorscaling = ["log";"scaled";"scaled";"scaled";"log";"log";"log"];
tbl = cat(1,splittbl{rbms});
colorlims    = [0      log(500)
                0.0224 0.4807
                0.0067 0.3569
                0.1250 0.3281
                log(0.0024) log(0.0772)
                log(0.0424) log(1.6010)
                -1 0];
               
%% Create heatmap
for mode = rbms
    %Create figure, figure appearance
    f{mode} = figure;
    clf
    sgtitle([figure_titles(criteria_idx),"RBM:"+str_mode(mode)]);
    f{mode}.Color = 'w';
    f{mode}.Position = [0 0 1920/8 1080];
    
    %create subplots for each sector angle
    for i = 2:6
        subplot(5,1,i-1)
        if criteria_idx == 1 || criteria_idx == 7
        h{mode,i} = heatmap(splittbl{mode},"R2",str_vari(i),'XDisplayData',XData,'YDisplayData',YData);
%         h{mode,i}.ColorLimits      = [0 log(500)];
%         elseif criteria_idx == 7
%         h{mode,i} = heatmap(splittbl{mode},"R2",str_vari(i),'XDisplayData',XData,'YDisplayData',YData);
        h{mode,i}.ColorLimits      = [0 1];
        else
        h{mode,i} = heatmap(splittbl{mode},"R2",str_vari(i),'ColorVariable',criterium        ,'ColorMethod',colormethod(criteria_idx),'XDisplayData',XData,'YDisplayData',YData);
        end
        
        %Labels
        h{mode,i}.Title            = heatmap_titles(criteria_idx);
        h{mode,i}.YLabel           = str_vari_text(i);
        h{mode,i}.MissingDataLabel = "-";
        h{mode,i}.FontName         = 'Times New Roman';
        
        %Ticks
        h{mode,i}.XDisplayLabels = XLabels;
        h{mode,i}.YDisplayLabels = YLabels;
        
        %Appearance
        if any(criteria_idx == 1:4)
            h{mode,i}.Colormap         = turbo;
        elseif any(criteria_idx == 1:6)
        	h{mode,i}.Colormap         = flipud(turbo);
        else
            h{mode,i}.Colormap         = gray;
        end

        h{mode,i}.MissingDataColor = 'k';
        h{mode,i}.GridVisible      = false;
        
        %Scaling
        h{mode,i}.ColorScaling     = colorscaling(criteria_idx);
        h{mode,i}.ColorLimits      = colorlims(criteria_idx,:);

    end
    %Save figures
    if nargin>3
    saveas(gcf, sprintf('Results/figures/%s_heatmap/%s_heatmap_%s_%s',runID,runID,criterium,str_mode(mode)), 'bmp')
    saveas(gcf, sprintf('Results/figures/%s_heatmap/%s_heatmap_%s_%s',runID,runID,criterium,str_mode(mode)), 'svg')
    saveas(gcf, sprintf('Results/figures/%s_heatmap/%s_heatmap_%s_%s',runID,runID,criterium,str_mode(mode)), 'pdf')
    else
        warning('no runID given, figures not saved')
    end
end
    