%% Gene association analysis
%% prepare predictor variable and response variable for PLS regression
data_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_data\';
load([data_dir,'prepro_gene\parcelExpression.mat']);
load([data_dir,'multimod_detecres\r1w1para\cross_sub_MV.mat']);
region_ind=find(~isnan(parcelExpression(:,2)));
response_var_file=cross_sub_MV(region_ind);
predictor_var_file=parcelExpression(region_ind,2:end);
MRIdata=response_var_file;
GENEdata=predictor_var_file;
geneindex=1:size(GENEdata,2);

%% PLS_calculation
Y=zscore(MRIdata);
dim=10;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(GENEdata,Y,dim,'CV',dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim); 

%align PLS components with desired direction%
R1 = corr([XS(:,1),XS(:,2),XS(:,3)],MRIdata);
if R1(1,1)<0
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0
    XS(:,2)=-1*XS(:,2);
end
if R1(3,1)<0
    XS(:,3)=-1*XS(:,3);
end

%calculate correlations of PLS components with MRI variables
corr_real = zeros(dim,1);
for cd=1:dim
    corr_real(cd) = corr(XS(:,cd),MRIdata);
end

%% validate the significance
load([data_dir,'SpinTest_result\finalres\permMV_res.mat']);
nperms=size(permMV_res,2);
PCTVARrand = zeros(nperms,dim);
perm_Rsquared = zeros(nperms,1);
for j=1:nperms
    disp(j);
    surro_Y=permMV_res(:,j);
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(GENEdata,zscore(surro_Y),dim);
     PCTVARrand(j,:)=PCTVARr(2,:);
     temp=cumsum(100*PCTVARr(2,1:dim));
     perm_Rsquared(j) = temp(dim); 
     % calculate correlations of PLS components with MRI variables
     for pd=1:dim
         [R_with_pls,~]=corr(XSr(:,pd),surro_Y);
         corr_permutate_PLS(pd,j) = R_with_pls;
     end
end
Covexp = zeros(1,dim);
for l=1:dim
    Covexp(l)=(length(find(PCTVARrand(:,l)>=PCTVAR(2,l)))+1)/(j+1);
end
p_cum =(length(find(perm_Rsquared>=Rsquared))+1)/(j+1);
[~,scovID]=sort(PCTVAR(2,:),'descend');
myStats=[PCTVAR(:,scovID); Covexp(:,scovID)];
csvwrite([data_dir,'genePLS_stats.csv'],myStats);

Corp = zeros(dim,1);
for d=1:dim
    Corp(d) =(length(find(corr_permutate_PLS(d,:)>=corr_real(d)))+1)/(j+1);
end
scorr_permutate_PLS=corr_permutate_PLS(scovID,:);
scorr_real=corr_real(scovID,:);
sCorp=Corp(scovID,:);

save([data_dir,'MV_corr_gene_permPLS.mat'],'scorr_permutate_PLS','scorr_real','sCorp');

%% Draw variance explanation
load([data_dir,'genePLS_stats.csv']);
py = plot(sort(genePLS_stats(2,:),'descend'),'.-');
py.Color = [0 0 0]/255;
py.MarkerSize = 10;
hold on

xlabel('PLS components');
ylabel('Variance explained (%)');
set(gca,'XLim',[0.5,10]);
set(gca,'YLim',[0.00,0.25],'YTick',0.00:0.05:0.25);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial','FontSize',10);
set(gca, 'TickLength', [0.02, 0.02],'TickDir', 'out');
box off
