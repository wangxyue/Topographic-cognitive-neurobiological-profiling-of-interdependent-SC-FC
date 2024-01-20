%% calculate corrected weight
data_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_data\';
load([data_dir,'prepro_gene\parcelExpression.mat']);
load([data_dir,'prepro_gene\probeInformation.mat'])
load([data_dir,'multimod_detecres\r1w1para\cross_sub_MV.mat']);
region_ind=find(~isnan(parcelExpression(:,2)));
response_var_file=cross_sub_MV(region_ind);
predictor_var_file=parcelExpression(region_ind,2:end);
MRIdata=response_var_file;
GENEdata=predictor_var_file;
geneindex=1:size(GENEdata,2);
gene_name=probeInformation.GeneSymbol;

%% PLS_calculation 
X=GENEdata;
Y=zscore(MRIdata);
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

load([data_dir,'genePLS_stats.csv']);
load([data_dir,'MV_corr_gene_permPLS.mat']);
sigdim=intersect(find(genePLS_stats(3,:)<0.05),find(Corp<0.05));
[R,p]=corr([XS(:,sigdim)],MRIdata);
if R(1,1)<0
    stats.W(:,sigdim)=-1*stats.W(:,sigdim);
    XS(:,sigdim)=-1*XS(:,sigdim);
end
[PLSw,x] = sort(stats.W(:,sigdim),'descend');
PLSids=gene_name(x);
PLSgeneindexp=geneindex(x);
PLS_ROIscores_360=zeros(360,1);
PLS_ROIscores_360(region_ind)=XS(:,sigdim);

save([data_dir,'genePLS_ROIscore.mat'],'PLS_ROIscores_360')

%% start bootstrap
disp(' Bootstrapping ')
%number of bootstrap iterations
bootnum=1000;
%define variables for storing the (ordered) weights from all bootstrap runs
boot_PLSweights=[];
for i=1:bootnum
    disp(i)
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled regions
    Yr=Y(myresample,:); % define Y for resampled regions
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
    
    temp=stats.W(:,sigdim);%extract PLS weights
    newW=temp(x); %order the newly obtained weights the same way as initial PLS 
    if corr(PLSw,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    boot_PLSweights=[boot_PLSweights,newW]; %store (ordered) weights from this bootstrap run    
end

% get standard deviation of weights from bootstrap runs
PLSsw=std(boot_PLSweights');
% get bootstrap weights
boot_ratio=PLSw./PLSsw';
boot_ratio(isnan(boot_ratio))=0;

%% order bootstrap weights and names of regions
[Zp_d,indp_d]=sort(boot_ratio,'descend');
PLSp_d=PLSids(indp_d);
geneindexp_d=PLSgeneindexp(indp_d);
fid1 = fopen(fullfile(data_dir,'PLS_geneWeights_descend.csv'),'w');
for i=1:length(gene_name)
  fprintf(fid1,'%s, %d, %f\n', PLSp_d{i},geneindexp_d(i), Zp_d(i));
end
fclose(fid1);

[Zp_a, indp_a]=sort(boot_ratio,'ascend');
PLSp_a=PLSids(indp_a);
geneindexp_a=PLSgeneindexp(indp_a);
fid2 = fopen(fullfile(data_dir,'PLS_geneWeights_ascend.csv'),'w');
for j=1:length(gene_name)
    fprintf(fid2,'%s, %d, %f\n', PLSp_a{j},geneindexp_a(j),Zp_a(j));
end
fclose(fid2);