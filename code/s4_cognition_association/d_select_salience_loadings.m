%% The reliability of variable loadings
% define work dir
data_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_data\';
load([data_dir,'hetepara_cogres.mat']);
load([data_dir,'HCPcogdata\Cognition_name.mat']);
load([data_dir,'hetepara_node.mat']);
sig_LC=intersect(find(hetepara_res.LC_pvals<0.05),find(hetepara_res.corr_LxLy_pvals<0.05));

%% select cognitive terms with salience weights
domains=Cognition_name;
beh_loadings=hetepara_res.LC_behav_loadings;
for nLC=1:size(beh_loadings,2)
    for beh=1:size(beh_loadings,1)
        L_CI=hetepara_res.boot_results.LC_behav_loadings_lB(beh,nLC);
        U_CI=hetepara_res.boot_results.LC_behav_loadings_uB(beh,nLC);
        if 0>=L_CI&&0<=U_CI
            beh_loadings(beh,nLC)=0;
        end
    end
end
for i=1:length(sig_LC)
    L=sig_LC(i);
    clear sig_beh_index sig_beh_loadings sig_beh_name ranked_sig_beh_loadings ranked_sig_ID
    sig_beh_index=find(beh_loadings(:,L)~=0);
    sig_beh_loadings=beh_loadings(sig_beh_index,1);
    sig_beh_name=Cognition_name(sig_beh_index);
    [ranked_sig_beh_loadings,ranked_sig_ID]=sort(sig_beh_loadings,'descend');
    final_sigbeh_loadings{1,i}=ranked_sig_beh_loadings;
    final_sigbeh_name{1,i}=sig_beh_name(ranked_sig_ID)';
end

save([data_dir,'hetepara_cogsigbeh.mat'],'final_sigbeh_name','final_sigbeh_loadings');

%% select brain regions with salience weights
MV_loadings=hetepara_res.LC_img_loadings;
for nLC=1:size(MV_loadings,2)
    for node=1:size(MV_loadings,2)
        L_CI=hetepara_res.boot_results.LC_img_loadings_lB(node,nLC);
        U_CI=hetepara_res.boot_results.LC_img_loadings_uB(node,nLC);
        if 0>=L_CI&&0<=U_CI
            MV_loadings(node,nLC)=0;
        end
    end
end
for i=1:length(sig_LC)
    L=sig_LC(i);
    clear sig_node_index
    sig_node_index=find(MV_loadings(:,L)~=0);
    final_MVloadings=zeros(360,1);
    final_MVloadings(hetepara_node,1)=MV_loadings(:,L);
    final_sigMV_loadings{1,i}=final_MVloadings;
end
save([data_dir,'hetepara_cogsigMV_loadings.mat'],'final_sigMV_loadings');
%surf_plot(final_MVloadings,'mycolor_bluered4.mat');
