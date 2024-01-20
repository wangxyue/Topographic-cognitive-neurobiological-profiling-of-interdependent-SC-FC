%% Regional multilayer modular variability
%% define work dir
data_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_data\';
load([data_dir,'sub1012_multilayer_net.mat']);
addpath(genpath(strcat(data_dir,'functions\GenLouvain2.1')));
intralayer_resolution = 1; 
interlayer_resolution = 1; 
layer_number=size(multilayer_net,2);
N_nodes=360;
run_times=100; 
%% Multilayer modularity detection & MV
for sub=1:size(multilayer_net,1) 
    for layer=1:layer_number
        AA{layer}=multilayer_net{sub,layer}; 
    end
    for i_time = 1:run_times
        [B,mm] = multiord(AA,intralayer_resolution,interlayer_resolution);
        PP = @(S) postprocess_ordinal_multilayer(S,layer_number); 
        [S,Q1,mod_iter_tmp] =iterated_genlouvain(B,10000,0,1,'moverandw',[],PP); 
        S = reshape(S,N_nodes,layer_number);
        Q_value_all(i_time,1) = Q1/mm;
        modularity_iterations_all(i_time,1) = mod_iter_tmp; 
        number_of_communities_all(i_time,1) = length(unique(S));
        modular_result_time{i_time,1} = S; 
        modular_result_transport_time{i_time,1}=modular_result_time{i_time,1}'; 
        modres = modular_result_transport_time{i_time,1};
        modvar = scaled_inclusivity_wei(modres); 
        MV_all(i_time,:) = modvar';
    end
    MV_mean = mean(MV_all);
    all_sub_MV(:,sub)=MV_mean';
    sub_multi_module{sub,1}=modular_result_time;
    sub_modularity{sub,1}=Q_value_all;
end
cross_sub_MV=mean(all_sub_MV,2);
% plot group-MV map
surf_plot(cross_sub_MV,'bocolorbar.mat')

save([data_dir,'multimod_detecres\r1w1para\all_sub_MV.mat'],'all_sub_MV');
save([data_dir,'multimod_detecres\r1w1para\sub_multi_module.mat'],'sub_multi_module');
save([data_dir,'multimod_detecres\r1w1para\sub_modularity.mat'],'sub_modularity');
save([data_dir,'multimod_detecres\r1w1para\cross_sub_MV.mat'],'cross_sub_MV');

%% Correlation with cortical gradient
load([data_dir,'multimod_detecres\r1w1para\cross_sub_MV.mat']);
load([data_dir,'gradient_glasser360.mat']);
[grad_r_real,~]=corr(cross_sub_MV,gradient_glasser360,'type','pearson');

% Compute the significance
load([data_dir,'SpinTest_result\finalres\Perm_360_MVResults.mat']);
for perm=1:size(Perm_360_Results,2)
    perm_MV=Perm_360_Results(:,perm);
    r=corr(perm_MV,gradient_glasser360,'type','pearson');
    grad_r_perm(perm)=r;
end
grad_sig_p=(length(find(grad_r_perm>=grad_r_real))+1)/(perm+1);
save([data_dir,'corr_grad.mat'],'grad_r_perm','grad_r_real','grad_sig_p');

% plot rvalue
h = histogram(grad_r_perm,60);
h.FaceColor = [0.3 0.3 0.3]; 
h.EdgeColor = [0.3 0.3 0.3];
set(gca, 'XTick', []);
set(gca, 'YTick', []);
box off
hold on
line([grad_r_real grad_r_real],[0 600], 'LineStyle', '--', 'Color', 'r');

%% Correlation with cortical expansion
load([data_dir,'multimod_detecres\r1w1para\cross_sub_MV.mat']);
load([data_dir,'CorticalExpansion_R_glasser360.mat']);
MV_R=cross_sub_MV(1:180,1);
[expan_r_real,~]=corr(MV_R,CorticalExpansion_R_glasser360,'type','pearson');

% Compute the significance
load([data_dir,'SpinTest_result\finalres\Perm_360_MVResults.mat']);
Perm_R_MV=Perm_360_Results(1:180,:);
for perm=1:size(Perm_R_MV,2)
    perm_MV=Perm_R_MV(:,perm);
    r=corr(perm_MV,CorticalExpansion_R_glasser360,'type','pearson');
    expan_r_perm(perm)=r;
end
expan_sig_p=(length(find(expan_r_perm>=expan_r_real))+1)/(perm+1);
save([data_dir,'corr_expan.mat'],'expan_r_perm','expan_r_real','expan_sig_p');

% plot rvalue
h = histogram(expan_r_perm,60);
h.FaceColor = [0.3 0.3 0.3]; 
h.EdgeColor = [0.3 0.3 0.3];
set(gca, 'XTick', []);
set(gca, 'YTick', []);
box off
hold on
line([expan_r_real expan_r_real],[0 700], 'LineStyle', '--', 'Color', 'r');

