%% neurotransmitter receptors and transporters distributions association analysis
% define work dir
data_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_data\';
load([data_dir,'multimod_detecres\r1w1para\cross_sub_MV.mat']);
load([data_dir,'SpinTest_result\finalres\Perm_360_MVResults.mat']);
neurotransmitter={'Acetylcholine','Cannabinoid','Dopamine','GABA','Glutamate','Histamine','Norepinephrine','Opioid','Serotonin'};

cormat=[19 5];
varTypes={'string','string','double','double','double'};
varNames={'neurotransmitter','receptors/transporters','real_r','sig_p','FDR_p'};
T=table('Size',cormat,'VariableTypes',varTypes,'VariableNames',varNames);
t=0;

for r=1:length(neurotransmitter)
    clear data_list data_num
    file_path=strcat(data_dir,'receptordata\nine_neurotransmitter_systems\',neurotransmitter{r},'\');
    data_list=dir(strcat(file_path,'*.mat'));
    data_num=length(data_list);
    for i=1:data_num
        clear rece_info r_real perm_MV r_perm sig_p
        rece_info=cell2mat(struct2cell(load(strcat(file_path,data_list(i).name))));
        [r_real,~]=corr(cross_sub_MV,rece_info,'type','pearson');
        
        %validation the significance
        for perm=1:size(Perm_360_Results,2)
            perm_MV=Perm_360_Results(:,perm);
            pr=corr(perm_MV,rece_info,'type','pearson');
            r_perm(perm)=pr;
        end
        sig_p=(length(find(r_perm>r_real))+1)/(perm+1);
        t=t+1;
        T(t,1:4)={neurotransmitter{r},data_list(i).name(1:end-4),r_real,sig_p};
        all_perm{t,1}=r_perm;
    end
end
%% FDR correction
ori_p=T.sig_p;
FDR_P=mafdr(ori_p,'BHFDR', true);
T.FDR_p=FDR_P;

[s_real_r,s_ID]=sort(T.real_r,'descend');
new_sortT=T(s_ID,:);

sigindex=find(FDR_P<0.05);
final_sig_T=T(sigindex,:);
final_perm=all_perm(sigindex);

writetable(new_sortT,[data_dir,'MV_cor_receptors.csv']);
save([data_dir,'MVcorrecep_resT.mat'],'new_sortT','final_sig_T','final_perm');

%% plot receptors maps
path=strcat(data_dir,'receptordata\nine_neurotransmitter_systems\');
load(strcat(path,'Acetylcholine\A4B2.mat'));
Z_A4B2=inormal(A4B2);
surf_plot(Z_A4B2,'bwocolorbar.mat')

load(strcat(path,'Cannabinoid\CB1.mat'));
Z_CB1=inormal(CB1);
surf_plot(Z_CB1,'bwocolorbar.mat')

load(strcat(path,'Opioid\MOR.mat'));
Z_MOR=inormal(MOR);
surf_plot(Z_MOR,'bwocolorbar.mat')

load(strcat(path,'Serotonin\R5HT4.mat'));
Z_R5HT4=inormal(R5HT4);
surf_plot(Z_R5HT4,'bwocolorbar.mat')
