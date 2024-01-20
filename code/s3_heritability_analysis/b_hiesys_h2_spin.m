%% Spintest
%% define work dir
data_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_data\';
figure_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_resfigs\';
NumPerm=10000; 
Perm_360_h2Results = zeros(360,NumPerm);
load('node_heritability.mat');
Perm_data = zscore(node_heritability);

for i=1:NumPerm
    disp(i)
    load(OutFile);
    Rand_Label = [LNewLabel;RNewLabel];
    for node=1:360
        index = find(GlasserLR==node);
        perm_roi = mode(Rand_Label(index));
        if perm_roi == 0
            Perm_360_h2Results(node,i)=0;
            continue;
        end
        Perm_360_h2Results(node,i) =Perm_data(perm_roi);
    end
end
save([data_dir,'Spintest\SpinTest_result\finalres\Perm_360_h2Results.mat'],'Perm_360_h2Results');


%% for 4 hier
load([data_dir,'SystemID_in_Glasser360\hier4_atlas.mat']);
null_models = zeros(NumPerm,4);
for i = 1:NumPerm
    tmp = Perm_360_h2Results(:,i);
    null_models(i,1) = mean(tmp(final_hier_360==1));
    null_models(i,2) = mean(tmp(final_hier_360==2));
    null_models(i,3) = mean(tmp(final_hier_360==3));
    null_models(i,4) = mean(tmp(final_hier_360==4));
end

for i = 1:4
    mean_hier4(i) = mean(Perm_data(final_hier_360==i));% mean value in 4 hier
end

for i = 1:4
    % mean R-square for real community ii
    x = mean_hier4(i);
    % mean R-square for permuted community ii
    mu = mean(squeeze(null_models(:,i)));
    % standard devitation of R-square for permuted community ii
    sigma = std(squeeze(null_models(:,i)));
    % z-score
    z_values(i) = (x - mu) / sigma;
end

for i=1:4
    data = null_models(:,i);
    if mean_hier4(i)>0
        mean_hier4_p(i) =  numel(find(data > mean_hier4(i))) / NumPerm;
    else
        mean_hier4_p(i) =  numel(find(data < mean_hier4(i))) / NumPerm;
    end
end

h2_data_hier4 = [z_values;mean_hier4_p];
save([data_dir,'index_spinres\h2_data_hier4.mat'],'h2_data_hier4');

%% plot sysh2
X = h2_data_hier4(1,:); 
hold on
color_matrix = [55,103,149;139,137,137;139,137,137;139,137,137]./255; 
for i = 1:4
    b = bar(i,X(i),'stacked','BarWidth',0.7);
    set(b(1),'facecolor',color_matrix(i,:))
end
box off
Xlabel = {'Pri', 'Uni', 'Hete', 'Para'};
set(gca,'XTick',1:length(X));
set(gca,'XTickLabel',Xlabel,'FontSize', 9);
xtickangle(45);
set(gca,'YLim',[-2,3],'YTick',-2:1:3);
ylabel('Z') 
set(gca,'FontSize',15,'Fontname', 'Arial');
text(1,X(1),'*','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',15,'FontName','Arial');
set(gca, 'TickDir', 'out','TickLength', [0.02, 0.025]);
set(gca, 'LineWidth', 1);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[0 0 5 7]);
print(gcf,[figure_dir,'sysh2_Z.tif'],'-dtiff','-r300')