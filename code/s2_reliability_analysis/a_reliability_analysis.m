%% test-retest reliability of multilayer modular variability
% define work dir
data_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_data\';
load([data_dir,'multimod_detecres\r1w1para\sub42_MV_test.mat']);
load([data_dir,'multimod_detecres\r1w1para\sub42_MV_retest.mat']);
modvar_test=sub_MV_test;
modvar_retest=sub_MV_retest;
sub_num=size(modvar_test,1);
N_nodes=360;
for i=1:sub_num
    for j=1:sub_num
        r=corr(modvar_test{i,1},modvar_retest{j,1},'type','pearson');
        MV_R(j,1)=r;
    end
    all_MV_R(:,i)=MV_R; 
end
real_intrasub_similarity=diag(all_MV_R);
T=all_MV_R;
T(logical(eye(size(T))))=0;
inter_sub=sum(T,1)/(sub_num-1);
real_intersub_similarity=inter_sub';
real_diff=mean(real_intrasub_similarity-real_intersub_similarity);
save([data_dir,'real_intrasub_similarity.mat'],'real_intrasub_similarity');
save([data_dir,'real_intersub_similarity.mat'],'real_intersub_similarity');

%% nonparametric permutation test 
nperms=10000;
for perm=1:nperms
    randind1=randperm(sub_num);
    for x=1:length(randind1)
        perm_MVtest{x,1}=modvar_test{randind1(x),1};
    end
    randind2=randperm(sub_num);
    for y=1:length(randind2)
        perm_MVretest{y,1}=modvar_retest{randind2(y),1};
    end
    for i=1:sub_num
        for j=1:sub_num
            perm_r=corr(perm_MVtest{i,1},perm_MVretest{j,1},'type','pearson');
            perm_MV_R(j,1)=perm_r;
        end
        all_perm_MV_R(:,i)=perm_MV_R;
    end
    perm_intra_sub=diag(all_perm_MV_R);
    permT=all_perm_MV_R;
    permT(logical(eye(size(permT))))=0;
    perm_inter_sub=sum(permT,1)/(sub_num-1);
    perm_inter_sub=perm_inter_sub';

    perm_intrasub_similarity(:,perm)=perm_intra_sub;
    perm_intersub_similarity(:,perm)=perm_inter_sub;
    
    all_perm_PR{perm,1}=all_perm_MV_R;
end
perm_diff=mean((perm_intrasub_similarity-perm_intersub_similarity),1);
diff_pvalue=length(find(perm_diff>=real_diff))/nperms;
save([data_dir,'reliability_res.mat'],'perm_diff','real_diff','diff_pvalue');

%% intra-class correlation(ICC) of multilayer modular variability
%poor: ICC<0.4; 
%moderate: 0.4<=ICC<0.6; 
%good: 0.6<=ICC<0.75; 
%excellent: ICC>=0.75
for node=1:N_nodes
    for sub=1:sub_num
        T_R_MV(sub,1)=modvar_test{sub,1}(node,1);
        T_R_MV(sub,2)=modvar_retest{sub,1}(node,1);
    end
    node_icc=Liao_icc(T_R_MV);
    ICC(node,1)=node_icc;
end
save([data_dir,'ICC.mat'],'ICC');

%% correlation between MV and ICC
load([data_dir,'multimod_detecres\r1w1para\cross_sub_MV.mat']);
load([data_dir,'Spintest\SpinTest_result\finalres\Perm_360_MVResults.mat']);
r_real=corr(ICC,cross_sub_MV,'type','pearson');
for perm=1:size(Perm_360_Results,2)
    perm_MV=Perm_360_Results(:,perm);
    r=corr(perm_MV,ICC,'type','pearson');
    r_perm(perm)=r;
end
r_sig_p=(length(find(r_perm>=r_real))+1)/(perm+1);
save([data_dir,'corMVICC_res.mat'],'r_perm','r_real','r_sig_p');
inor_MV=inormal(cross_sub_MV);
inor_ICC=inormal(ICC);

%% plot rvalue
h = histogram(r_perm,60);
h.FaceColor = [0.3 0.3 0.3]; 
h.EdgeColor = [0.3 0.3 0.3];
set(gca, 'XTick', []);
set(gca, 'YTick', []);
box off
hold on
line([r_real r_real],[0 800], 'LineStyle', '--', 'Color', 'r','LineWidth', 2);