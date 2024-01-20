%% Cognition association analysis
% define work dir
data_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_data\';
load([data_dir, 'HCPcogdata\ori_Cognition.mat']);
load([data_dir, 'multimod_detecres\r1w1para\priuni_MV.mat']);
X0=priuni_MV;
Y0=ori_Cognition;
[nSubj,node] = size(X0); 
ncog_terms = size(Y0,2);
% Normalize input data matrices
X=zscore(X0);
Y=zscore(Y0);

%% Cross-covariance matrix
R=Y'*X;
% Singular value decomposition
[U,S,V] = svd(R,'econ');
% Number of latent components (LCs)
nLC = min(size(S)); 
% ICA convention: turn LCs such that max is positive
for iLC = 1:nLC
    [~,maxID] = max(abs(V(:,iLC)));
    if sign(V(maxID,iLC))<0
        V(:,iLC) = -V(:,iLC);
        U(:,iLC) = -U(:,iLC);
    end
end
% Amount of covariance explained by each LC
explCovLC = (diag(S).^2) / sum(diag(S.^2));
%% Compute PLS scores & loadings
Lx = X * V;
Ly = Y * U;
corr_Lx_X = corr(Lx,X)';
corr_Ly_Y = corr(Ly,Y)';

for num_LC=1:nLC
    cor_real_LxLy(:,num_LC)=corr(Lx(:,num_LC),Ly(:,num_LC));
end

%% Permutation testing for LC significance
% Set up random number generator
rng(1)
disp('... Permutations ...')
nPerms=10000;
for iP = 1:nPerms
    disp(iP);
    Xp = X;
    Yp = zeros(size(Y));
    perm_order = randperm(size(Y,1));
    Yp = Y(perm_order,:);
    % Generate cross-covariance matrix between Xp and permuted Y
    Rp =Yp'*Xp;
    [Up,Sp,Vp] = svd(Rp,'econ');
    % Compute PLS scores & loadings
    Lxp = Xp * Vp;
    Lyp = Yp * Up;
    PLS_xscore{iP,1}=Lxp;
    PLS_yscore{iP,1}=Lyp;
    corr_Lxp_Xp = corr(Lxp,Xp)';
    corr_Lyp_Yp = corr(Lyp,Yp)';
    PLS_LxXload{iP,1}=corr_Lxp_Xp;
    PLS_LyYload{iP,1}=corr_Lyp_Yp;
    % Procrustas transform (correction for axis rotation/reflection)
    rotatemat = rri_bootprocrust(U, Up);
    Up = Up * Sp * rotatemat;
    Sp = sqrt(sum(Up.^2));
    % Keep singular values for sample distribution of singular values
    Sp_vect(:,iP) = Sp'; 
end
%% Compute the significance from the permutation null distribution
S_mat = repmat(diag(S),1,nPerms);
sp = sum(Sp_vect >= S_mat,2);
LC_pvals = (sp + 1) ./ (nPerms + 1);
signif_LC = find(LC_pvals<0.05); % index of significant LCs
nSignifLC = size(signif_LC,1); % number of significant LCs
for iLC = 1:nSignifLC
    this_lc = signif_LC(iLC);
    disp(['LC' num2str(this_lc) '-p=' num2str(LC_pvals(this_lc),'%0.3f') ]);
end

% calculate the corr(Lxp,Lyp) from the permutation null distribution
for run=1:nPerms
    for num_LC=1:nLC
        perm_Lx=PLS_xscore{run,1}(:,num_LC);
        perm_Ly=PLS_yscore{run,1}(:,num_LC);
        cor_perm_LxLy(run,num_LC)=corr(perm_Lx,perm_Ly);
    end
end
% calculate the significance of corr(Lx,Ly) for each LC
for num_LC=1:nLC
    perm_value=cor_perm_LxLy(:,num_LC);
    real_value=cor_real_LxLy(num_LC);
    c=length(find(perm_value>=real_value));
    corr_LxLy_pvals(num_LC) = (c + 1) ./ (nPerms + 1);
end

%% Bootstrapping to test stability of PLS loadings
disp('... Bootstrapping ...')
% Set up random number generator
rng(1);
nBootstraps=1000;
% Get bootstrap subject sampling
all_boot_orders = nan(nSubj,nBootstraps);
[boot_order_tmp,~] = rri_boot_order(nSubj,1,nBootstraps);
all_boot_orders= boot_order_tmp;
% Run PLS in each bootstrap sample
for iB = 1:nBootstraps
    Xb = X0(all_boot_orders(:,iB),:);
    Xb = zscore(Xb);
    Yb = Y0(all_boot_orders(:,iB),:);
    Yb = zscore(Yb);
    Rb =Yb'*Xb;
    [Ub,Sb,Vb] = svd(Rb,'econ');
    % Procrustas transform (correction for axis rotation/reflection)
    % Computed on U only
    rotatemat_full = rri_bootprocrust(U, Ub);
    % Rotate and re-scale Ub and Vb
    Vb = Vb * Sb * rotatemat_full;
    Ub = Ub * Sb * rotatemat_full;
    Vb = Vb./repmat(diag(S)',size(Vb,1),1);
    Ub = Ub./repmat(diag(S)',size(Ub,1),1);
    % Vectors with all bootstrap samples -> needed for percentile computation
    Ub_vect(:,:,iB) = Ub;
    Vb_vect(:,:,iB) = Vb;
    % Compute bootstrapping PLS scores
    Lxb=Xb * Vb;
    Lyb=Yb * Ub;
    corr_Lxb_Xb=corr(Lxb,Xb)';
    corr_Lyb_Yb=corr(Lyb,Yb)';
    boot_results.Lxb(:,:,iB) = Lxb;
    boot_results.Lyb(:,:,iB) = Lyb;
    boot_results.LC_img_loadings_boot(:,:,iB) = corr_Lxb_Xb;
    boot_results.LC_behav_loadings_boot(:,:,iB) = corr_Lyb_Yb;
end
% Compute bootstrapping statistics
boot_stats = myPLS_bootstrap_stats(Ub_vect,Vb_vect,boot_results);
% Save all the statistics fields in the boot_results 
fN = fieldnames(boot_stats);
for iF = 1:length(fN)
    boot_results.(fN{iF}) = boot_stats.(fN{iF});
end
boot_results.Ub_vect = Ub_vect;
boot_results.Vb_vect = Vb_vect;

%% Save all result variables in struct
priuni_res.X0 = X0;
priuni_res.Y0 = Y0;
priuni_res.X = X;
priuni_res.Y = Y;

priuni_res.R = R;
priuni_res.U = U;
priuni_res.S = S;
priuni_res.V = V;

priuni_res.explCovLC = explCovLC;
priuni_res.LC_pvals = LC_pvals;

priuni_res.Lx = Lx;
priuni_res.Ly = Ly;

priuni_res.cor_perm_LxLy=cor_perm_LxLy;
priuni_res.cor_real_LxLy=cor_real_LxLy;
priuni_res.corr_LxLy_pvals=corr_LxLy_pvals;

priuni_res.LC_img_loadings = corr_Lx_X; 
priuni_res.LC_behav_loadings = corr_Ly_Y;

priuni_res.Sp_vect = Sp_vect;
priuni_res.LC_pvals = LC_pvals;
priuni_res.boot_results = boot_results;

save([data_dir,'priuni_cogres.mat'],'priuni_res');

%% Draw covariance explanation
py = plot(sort(explCovLC(1:10),'descend'),'.-');
py.Color = [0 0 0]/255;
py.LineWidth = 1;
py.MarkerSize = 15;
hold on

xlabel('Latent variable');
ylabel('Covariance explained (%)');
set(gca,'XLim',[0.5,10]);
set(gca,'YLim',[0.0,0.6],'YTick',0.0:0.1:0.6);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial','FontSize',10);
set(gca, 'TickLength', [0.02, 0.02],'TickDir', 'out');
box off

