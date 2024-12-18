function ACEfit_Results(ACEfit_Par,varargin)
% FORMAT ACEfit_Results(ACEfit_Par[,FWEalpha,FDRalpha])
%_______________________________________________________________________
% Version: http://github.com/nicholst/APACE/tree/34ea883
%          2018-09-21 09:59:07 +0100

if nargin>=2
    FWEalpha = varargin{1};
else
    FWEalpha = 0.05;
end
if nargin>=3
    FDRalpha = varargin{2};
else
    FDRalpha = 0.05;
end

%
% Plotting Function
%

load(fullfile(ACEfit_Par.ResDir,'ACEfit_Perm'));
N             = ACEfit_Par.nPerm+1;
ctl_val_index = N-floor(N*FWEalpha); % floor(N*(1-FWEalpha))+1

f = 0;

%%% (1) MEAN
f = f+1;
SetFig(f);
[sMEAN,oMEAN] = sort(mean_ACE);
ctl_val_MEAN  = sMEAN(ctl_val_index);
[f1,x1]       = hist(sMEAN,100);
n_MEAN        = min(find(sMEAN(:)==mean_ACE(N)));
bar(x1,f1/N);

% Calculate the permutation-based p-value
p_MEAN = (N-n_MEAN+1)/N;

title(sprintf('H0 dist of Mean of h^2, P-value=%.3f',p_MEAN));
xlabel('mean of h^2');
yLimits = get(gca,'YLim');
line([mean_ACE(N) mean_ACE(N)],[0 yLimits(2)],'Marker','.','Color','green');
% Plot the critical threshold (red)
line([ctl_val_MEAN ctl_val_MEAN],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
text(mean_ACE(N),0.8*yLimits(2),num2str(mean_ACE(N)))
text(ctl_val_MEAN,0.9*yLimits(2),num2str(ctl_val_MEAN))
print('-dpdf',fullfile(ACEfit_Par.ResDir,'H0dist_mean.pdf'));

%%% (2) wh2
f = f+1;
SetFig(f);
[swh2,owh2] = sort(wh2_ACE);
ctl_val_wh2 = swh2(ctl_val_index);
[f1,x1]     = hist(swh2,100);
n_wh2       = min(find(swh2(:)==wh2_ACE(N)));
bar(x1,f1/N);

% Calculate the permutation-based p-value
p_wh2 = (N-n_wh2+1)/N;

title(sprintf('H0 dist of Weighted Mean of h^2, P-value=%.3f',p_wh2));
xlabel('weighted mean of h^2');
yLimits = get(gca,'YLim');
line([wh2_ACE(N) wh2_ACE(N)],[0 yLimits(2)],'Marker','.','Color','green');
% Plot the critical threshold (red)
line([ctl_val_wh2 ctl_val_wh2],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
text(wh2_ACE(N),0.8*yLimits(2),num2str(wh2_ACE(N)))
text(ctl_val_wh2,0.9*yLimits(2),num2str(ctl_val_wh2))
print('-dpdf',fullfile(ACEfit_Par.ResDir,'H0dist_wh2.pdf'));


%%% (3) MEDIAN
f = f+1;
SetFig(f);
[sMEDIAN,oMEDIAN] = sort(med_ACE);
ctl_val_MEDIAN    = sMEDIAN(ctl_val_index);
[f1,x1]           = hist(sMEDIAN,100);
n_MEDIAN          = min(find(sMEDIAN(:)==med_ACE(N)));
bar(x1,f1/N);

% Calculate the permutation-based p-value
p_MEDIAN = (N-n_MEDIAN+1)/N;

title(sprintf('H0 dist of Median (Q2) of h^2, P-value=%.3f',p_MEDIAN));
xlabel('Q2 of h^2');
yLimits = get(gca,'YLim');
line([med_ACE(N) med_ACE(N)],[0 yLimits(2)],'Marker','.','Color','green');
% Plot the critical threshold (red)
line([ctl_val_MEDIAN ctl_val_MEDIAN],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
text(med_ACE(N),0.8*yLimits(2),num2str(med_ACE(N)))
text(ctl_val_MEDIAN,0.9*yLimits(2),num2str(ctl_val_MEDIAN))
print('-dpdf',fullfile(ACEfit_Par.ResDir,'H0dist_median.pdf'));


%%% (4) Q3
f = f+1;
SetFig(f);
[sQ3,oQ3]  = sort(q3_ACE);
ctl_val_Q3 = sQ3(ctl_val_index);
[f1,x1]    = hist(sQ3,100);
n_Q3       = min(find(sQ3(:)==q3_ACE(N)));
bar(x1,f1/N);

% Calculate the permutation-based p-value
p_Q3 = (N-n_Q3+1)/N;

title(sprintf('H0 dist of Third Quartile (Q3) of h^2, P-value=%.3f',p_Q3));
xlabel('Q3 of h^2');
yLimits = get(gca,'YLim');
line([q3_ACE(N) q3_ACE(N)],[0 yLimits(2)],'Marker','.','Color','green');
% Plot the critical threshold (red)
line([ctl_val_Q3 ctl_val_Q3],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
text(q3_ACE(N),0.8*yLimits(2),num2str(q3_ACE(N)))
text(ctl_val_Q3,0.9*yLimits(2),num2str(ctl_val_Q3))
print('-dpdf',fullfile(ACEfit_Par.ResDir,'H0dist_q3.pdf'));


%%% (5) mean(>=median)
if any(isnan(mGmed_ACE))
    
    p_mGTmedian = NaN;
    
else
    
    f = f+1;
    SetFig(f);
    [smGTmedian,omGTmedian] = sort(mGmed_ACE);
    ctl_val_mGTmedian       = smGTmedian(ctl_val_index);
    [f1,x1]                 = hist(smGTmedian,100);
    n_mGTmedian             = min(find(smGTmedian(:)==mGmed_ACE(N)));
    bar(x1,f1/N);
    
    % Calculate the permutation-based p-value
    p_mGTmedian = (N-n_mGTmedian+1)/N;
    
    title(sprintf('H0 dist of Mean of h^2 >= Q2(h^2), P-value=%.3f',p_mGTmedian));
    xlabel('mean of h^2 > Q2(h^2)');
    yLimits = get(gca,'YLim');
    line([mGmed_ACE(N) mGmed_ACE(N)],[0 yLimits(2)],'Marker','.','Color','green');
    % Plot the critical threshold (red)
    line([ctl_val_mGTmedian ctl_val_mGTmedian],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
    text(mGmed_ACE(N),0.8*yLimits(2),num2str(mGmed_ACE(N)))
    text(ctl_val_mGTmedian,0.9*yLimits(2),num2str(ctl_val_mGTmedian))
    print('-dpdf',fullfile(ACEfit_Par.ResDir,'H0dist_mGTmedian.pdf'));
    
end


%%% (6) mean(>=q3)
if any(isnan(mGq3_ACE))
    
    p_mGTq3 = NaN;
    
else
    
    f = f+1;
    SetFig(f);
    [smGTq3,omGTq3] = sort(mGq3_ACE);
    ctl_val_mGTq3   = smGTq3(ctl_val_index);
    [f1,x1]         = hist(smGTq3,100);
    n_mGTq3         = min(find(smGTq3(:)==mGq3_ACE(N)));
    bar(x1,f1/N);
    
    % Calculate the permutation-based p-value
    p_mGTq3 = (N-n_mGTq3+1)/N;
    
    title(sprintf('H0 dist of Mean of h^2 >= Q3(h^2), P-value=%.3f',p_mGTq3));
    xlabel('mean of h^2 > Q3(h^2)');
    yLimits = get(gca,'YLim');
    line([mGq3_ACE(N) mGq3_ACE(N)],[0 yLimits(2)],'Marker','.','Color','green');
    % Plot the critical threshold (red)
    line([ctl_val_mGTq3 ctl_val_mGTq3],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
    text(mGq3_ACE(N),0.8*yLimits(2),num2str(mGq3_ACE(N)))
    text(ctl_val_mGTq3,0.9*yLimits(2),num2str(ctl_val_mGTq3))
    print('-dpdf',fullfile(ACEfit_Par.ResDir,'H0dist_mGTq3.pdf'));
    
end


Pvals_h2 = [p_MEAN p_wh2 p_MEDIAN p_Q3 p_mGTmedian p_mGTq3]';

save(fullfile(ACEfit_Par.ResDir,'Pvals_h2'),'Pvals_h2');


if ~ACEfit_Par.NoImg
    
    %%% (7) Maximum Test Statistic T
    f = f+1;
    SetFig(f);
    [sT,oT]   = sort(max_T_ACE);
    ctl_val_T = sT(ctl_val_index);
    [f3,x3]   = hist(sT,100);
    n_T       = min(find(sT(:)==max_T_ACE(N)));
    bar(x3,f3/N);
    
    % Calculate the permutation-based p-value
    p_T = (N-n_T+1)/N;
    
    title(sprintf('H0 dist of Maximum LRT Statistic (T), FWE P-value=%.3f',p_T));
    xlabel('T');
    ylabel('f(T)');
    yLimits = get(gca,'YLim');
    line([max_T_ACE(N) max_T_ACE(N)],[0 yLimits(2)],'Marker','.','Color','green');
    % Plot the critical threshold (red)
    line([ctl_val_T ctl_val_T],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
    text(max_T_ACE(N),0.8*yLimits(2),num2str(max_T_ACE(N)))
    text(ctl_val_T,0.9*yLimits(2),num2str(ctl_val_T))
    print('-dpdf',fullfile(ACEfit_Par.ResDir,'H0dist_test_statistic.pdf'));
    
    fprintf('The %.2f level critical threshold for the maximum LRT statistic is %8.3f. \n \n', FWEalpha, ctl_val_T);
    
    Vs = ACEfit_Par.Vs;
    
    if length(Vs.Dim)==1
        Dim = [Vs.Dim 1];
    else
        Dim = Vs.Dim;
    end
    
    %
    % Voxel-wise FDR correction
    %
    
    fdrPval = unPval_ACE(ACEfit_Par.I_data');
    nMsk = numel(fdrPval);
    
    % Compute FDR corrected p-values
    [sP,iP]  = sort(fdrPval);
    cV       = 1;
    Qs       = min(sP./(1:nMsk)'*nMsk*cV,1);
    sP(end)  = Qs(end);
    for iMsk = nMsk-1:-1:1
        sP(iMsk) = min(Qs(iMsk),sP(iMsk+1));
    end
    
    FDR_min  = min(Qs);
    fdrPval(iP) = sP;
    
    %%% (8) FDR plot
    f = f+1;
    SetFig(f);
    plot((1:nMsk)'/nMsk,sort(unPval_ACE(ACEfit_Par.I_data')),(1:nMsk)'/nMsk,sort(fdrPval));
    legend({'Uncorr P','FDR-corr P'},'Location','NorthWest')
    xlabel('Expected Uncorrected Ordered P');
    ylabel('Ordered P');
    title('Element-wise h^2 P-values: FDR plot');
    axis equal;axis tight;axis([0 1 0 1]);
    abline(0,1,'LineStyle','-','color',[.5 .5 .5]);
    abline(0,FDRalpha,'LineStyle',':','color','red')
    print('-dpdf',fullfile(ACEfit_Par.ResDir,'h2_FDRplot.pdf'));
    
    FDR_Pval = ones(Dim);
    FDR_Pval(ACEfit_Par.I_data') = fdrPval;
    
    % % Apply the threshold FDRalpha
    % FDR_Pval(FDR_Pval>FDRalpha) = 1;
    
    % -log10(FDR corrected p-value) for voxels
    FDR_Pval = -log10(FDR_Pval);
    
    fprintf('Minimum element-wise P_FDR is %.4f. \n \n', FDR_min);
    
    %
    % Voxel-wise FWE correction
    %
    
    % -log10(uncorrected P-value) for T
    unPval_ACE   = -log10(unPval_ACE);
    
    Tstat        = reshape(ACEfit_Par.Stats,Dim);
    corrPval_ACE = zeros(Dim);
    
    for i=1:numel(Tstat)
        if Tstat(i)>0
            % -log10(FWE corrected P-value) for T
            corrPval_ACE(i) = -log10(mean(Tstat(i)<=max_T_ACE));
        end
    end
    
    fprintf('Minimum element-wise P_FWE is %.4f. \n \n', 10.^(-max(corrPval_ACE(:))));

    
    %
    % Write out the output images
    %
    switch upper(ACEfit_Par.Model)
        case 'ACE'
            WriteData(unPval_ACE,   Vs, 'ACE_A_LRT_vox_P',    ACEfit_Par.ResDir);
            WriteData(corrPval_ACE, Vs, 'ACE_A_LRT_vox_FWEP', ACEfit_Par.ResDir);
            WriteData(FDR_Pval,     Vs, 'ACE_A_LRT_vox_FDRP', ACEfit_Par.ResDir);
        case 'AE'
            WriteData(unPval_ACE,   Vs, 'AE_A_LRT_vox_P',    ACEfit_Par.ResDir);
            WriteData(corrPval_ACE, Vs, 'AE_A_LRT_vox_FWEP', ACEfit_Par.ResDir);
            WriteData(FDR_Pval,     Vs, 'AE_A_LRT_vox_FDRP', ACEfit_Par.ResDir);
    end
    
    if ACEfit_Par.Vs.ClustInf
        
        %%% (9) Maximum Suprathreshold Cluster Size K
        f = f+1;
        SetFig(f);
        [sK,oK]   = sort(max_K_ACE);
        ctl_val_K = sK(ctl_val_index);
        [f1,x1]   = hist(sK,100);
        n_K       = min(find(sK(:)==max_K_ACE(N)));
        bar(x1,f1/N);
        
        % Calculate the permutation-based p-value
        p_K = (N-n_K+1)/N;
        
        title(sprintf('H0 dist of Maximum Suprathreshold Cluster Size (K), FWE P-value=%.3f',p_K));
        xlabel('K');
        ylabel('f(K)');
        yLimits = get(gca,'YLim');
        line([max_K_ACE(N) max_K_ACE(N)],[0 yLimits(2)],'Marker','.','Color','green');
        % Plot the critical threshold (red)
        line([ctl_val_K ctl_val_K],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
        text(max_K_ACE(N),0.8*yLimits(2),num2str(max_K_ACE(N)))
        text(ctl_val_K,0.9*yLimits(2),num2str(ctl_val_K))
        print('-dpdf',fullfile(ACEfit_Par.ResDir,'H0dist_cluster_size.pdf'));
        
        
        %%% (10) Maximum Suprathreshold Cluster Mass M
        f = f+1;
        SetFig(f);
        [sM,oM]   = sort(max_M_ACE);
        ctl_val_M = sM(ctl_val_index);
        [f2,x2]   = hist(sM,100);
        n_M       = min(find(sM(:)==max_M_ACE(N)));
        bar(x2,f2/N);
        
        % Calculate the permutation-based p-value
        p_M = (N-n_M+1)/N;
        
        title(sprintf('H0 dist of Maximum Suprathreshold Cluster Mass (M), FWE P-value = %.3f',p_M));
        xlabel('M');
        ylabel('f(M)');
        yLimits = get(gca,'YLim');
        line([max_M_ACE(N) max_M_ACE(N)],[0 yLimits(2)],'Marker','.','Color','green');
        % Plot the critical threshold (red)
        line([ctl_val_M ctl_val_M],[0 yLimits(2)],'Marker','.','LineStyle','-.','Color','red');
        text(max_M_ACE(N),0.8*yLimits(2),num2str(max_M_ACE(N)))
        text(ctl_val_M,0.9*yLimits(2),num2str(ctl_val_M))
        print('-dpdf',fullfile(ACEfit_Par.ResDir,'H0dist_cluster_mass.pdf'));
        
        fprintf('The %.2f level critical threshold for the maximum suprathreshold cluster size is %d. \n',    FWEalpha, ctl_val_K);
        fprintf('The %.2f level critical threshold for the maximum suprathreshold cluster mass is %8.3f. \n', FWEalpha, ctl_val_M);
        
        %
        % FWE correction
        %
        
        [ClSz,ClMass,corrPClSz,corrPClMass] = deal(zeros(Dim));
        
        % Cluster-forming threshold
        Tclus = spm_invXcdf(1-2*ACEfit_Par.alpha_CFT,1);
        % Apply Tclus threshold
        Tstat(Tstat(:)<Tclus) = 0;
        
        % Cluster-forming procedure
        [L,NUM] = spm_bwlabel(Tstat,18);
        
        if NUM>0
            cluster_size = zeros(NUM,1);
            cluster_mass = zeros(NUM,1);
            for i=1:NUM
                cluster_size(i) = sum(L(:)==i); % length(find(L(:)==i));
                cluster_mass(i) = sum(Tstat(L(:)==i));
                
                ClSz(L(:)==i)   = cluster_size(i);
                ClMass(L(:)==i) = cluster_mass(i);
                
                % -log10(FWE corrected P-value) for K
                corrPClSz(L(:)==i)   = -log10(mean(cluster_size(i)<=max_K_ACE));
                % -log10(FWE corrected P-value) for M
                corrPClMass(L(:)==i) = -log10(mean(cluster_mass(i)<=max_M_ACE));
            end
        end
        
        %
        % Write out the output images
        %
        switch upper(ACEfit_Par.Model)            
            case 'ACE'
                WriteData(ClSz,        Vs, 'ACE_A_LRT_clus',      ACEfit_Par.ResDir);
                WriteData(ClMass,      Vs, 'ACE_A_LRT_mass',      ACEfit_Par.ResDir);
                WriteData(corrPClSz,   Vs, 'ACE_A_LRT_clus_FWEP', ACEfit_Par.ResDir);
                WriteData(corrPClMass, Vs, 'ACE_A_LRT_mass_FWEP', ACEfit_Par.ResDir);
            case 'AE'
                WriteData(ClSz,        Vs, 'AE_A_LRT_clus',      ACEfit_Par.ResDir);
                WriteData(ClMass,      Vs, 'AE_A_LRT_mass',      ACEfit_Par.ResDir);
                WriteData(corrPClSz,   Vs, 'AE_A_LRT_clus_FWEP', ACEfit_Par.ResDir);
                WriteData(corrPClMass, Vs, 'AE_A_LRT_mass_FWEP', ACEfit_Par.ResDir);
        end
        
        save(fullfile(ACEfit_Par.ResDir,'Pvals_Max_h2'),'p_T','p_K','p_M');
        
    else
        
        save(fullfile(ACEfit_Par.ResDir,'Pvals_Max_h2'),'p_T');
        
    end
    
end

return

function SetFig(f)

figure(f)
set(gcf,'PaperPosition',[0 0 10 6])
set(gcf,'PaperSize',[10 6])

return
