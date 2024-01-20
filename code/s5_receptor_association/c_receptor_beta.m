%% plot receptor's beta
data_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_data\';
figure_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_resfigs\';
load([data_dir,'Elastic_beta.mat']);

[ns_beta,ind]=sort(final_beta,'descend');
ns_name=final_name(ind);
Bvalue=ns_beta';

figure(1)
RC = radarChart(Bvalue ,'Type','Patch');
RC.RLim = [-1,1];                         
RC.RTick = [-1,0.2,0.6,1.0];                   
RC.PropName = {'MOR','5-HT_4','\alpha_4\beta_2','5-HT_{1B}','DAT','CB_1','M_1','NET','5-HTT','H_3','GABA_A'};
RC.ClassName = {'\beta'};
RC = RC.draw(); 
RC.legend();                                 
colorList=[47 91 126;
          231 098 084;
          184 168 207;
          231 188 198;
          253 207 158;
          239 164 132;
          182 118 108]./255;
for n=1:RC.ClassNum
    RC.setPatchN(n,'FaceColor',colorList(n,:),'EdgeColor',colorList(n,:))
end
RC.setThetaTick('LineWidth',1,'Color',[211/255,211/255,211/255]);         
RC.setRTick('LineWidth',1,'Color',[030/255,070/255,110/255]);             
RC.setPropLabel('FontSize',14,'FontName','Arial','Color',[0,0,0])              
RC.setRLabel('FontSize',13,'FontName','Arial','Color',[.8,0,0])                    
% RC.setBkg('FaceColor',[0.8,0.8,0.8])                
% RC.setRLabel('Color','none')  

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 6 6]);
print(gcf,[figure_dir,'Elas_betavalue.tif'],'-dtiff','-r600')