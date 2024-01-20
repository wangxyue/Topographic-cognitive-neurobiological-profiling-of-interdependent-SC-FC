%% Elastic net regression model
% define work dir
data_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_data\';
figure_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_resfigs\';
filename = [data_dir, 'XYdata.csv'];
DATA=readtable(filename);
XYdata_matrix = table2array(DATA);
receptor_data=XYdata_matrix(:,2:end-1);
MV_data=XYdata_matrix(:,end);

% perform normalization
[m,c] = size(receptor_data); 
X = zeros(m,c); 
% normalize receptor data between -1 and 1
x=receptor_data;
for n=1:c
    X(:,n) = (x(:,n)-nanmean(x(:,n)))/(nanmax(x(:,n))-nanmin(x(:,n)));
end
% normalize MV data between -1 and 1
Y_tmp = MV_data; 
Y = (Y_tmp-nanmean(Y_tmp))/(nanmax(Y_tmp)-min(Y_tmp));

%% select optimal lambda
lambda = logspace(-5, 2, 100);
% alpha: 1 = lasso; 0 = ridge regression; in-between 0 and 1 = elastic net
alpha = 0.5; 
% elastic net regression with cross-validation of k=10
cv = cvpartition(size(X,1),'KFold',10);
[B,FitInfo] = lasso(X,Y,'Alpha',alpha,'NumLambda',length(lambda),'Lambda',lambda,'CV',cv); 

best_lambda=FitInfo.LambdaMinMSE;
Min_MSE=FitInfo.MSE(FitInfo.IndexMinMSE);
% beta values at minumum MSE based on cross-validation
B_minMSE = B(:,FitInfo.IndexMinMSE); 
real_beta=B_minMSE;
real_intercept=FitInfo.Intercept(FitInfo.IndexMinMSE);

nz_index=find(real_beta~=0);
load([data_dir,'receptordata\neurotrans_name.mat']);
final_name=neurotrans_name(nz_index,:);
final_beta=real_beta(nz_index);
save([data_dir,'Elastic_beta.mat'],'final_beta','final_name');

%% cross-validated Elastic net regression
lassoPlot(B,FitInfo, 'PlotType','CV' );
legend('off'); 
box off;
title('Cross-validated mean square of Elastic net fit');
h = gca;
set(h, 'TickDir', 'out','TickLength',[0.02 0.02]);
set(h, 'LineWidth', 1);
xlabel('\lambda', 'Interpreter', 'tex');
ylabel('Mean square error');
set(h,'YLim',[0.035,0.07],'YTick',0.035:0.005:0.07);

allLines = findobj(gca, 'type', 'line');
for i = 1:length(allLines)
    xdata = get(allLines(i), 'XData');
    if length(xdata) == 1 && abs(xdata - FitInfo.LambdaMinMSE) > eps
        delete(allLines(i));
    end
    if length(xdata) == 2 && abs(xdata(2) - FitInfo.LambdaMinMSE) > eps
        delete(allLines(i));
    end
end
lineX = [FitInfo.LambdaMinMSE FitInfo.LambdaMinMSE];
lineY = [h.YLim(1) h.YLim(2)];
line(lineX, lineY, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);
allMarkers = findobj(gca, 'Marker', 'o');
for i = 1:length(allMarkers)
    xdata = get(allMarkers(i), 'XData');
    if abs(xdata - FitInfo.LambdaMinMSE) < eps
        set(allMarkers(i), 'MarkerEdgeColor', 'k', 'LineStyle', '--'); 
    end
end

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 9 7]);
print(gcf,[figure_dir,'Elassel_optlambda.tif'],'-dtiff','-r600')

%% calculate model R2 (empirical)
yped=X(:,nz_index)*final_beta+real_intercept;
SS_Residual = sum((Y - yped).^ 2);
SS_Total = sum((Y - mean(Y)).^2);
real_model_Rsquare = 1 - (SS_Residual/ SS_Total);
real_rvalue=corr(Y,yped);
save([data_dir,'Elasreal_model_res.mat'],'real_model_Rsquare','real_rvalue');

%% validate the significance
load([data_dir,'SpinTest_result\finalres\Perm_360_MVResults.mat']);
spins=size(Perm_360_Results,2);
for j=1:spins
    disp(j);
    spin_Y=Perm_360_Results(:,j);
    rand_Y = (spin_Y-nanmean(spin_Y))/(nanmax(spin_Y)-min(spin_Y));
    [rand_B,rand_FitInfo] = lasso(X,rand_Y,'Alpha',alpha,'Lambda',best_lambda);

    random_y_pred =X(:,nz_index)*rand_B(nz_index)+rand_FitInfo.Intercept;
    rmSS_Residual = sum((rand_Y - random_y_pred).^ 2);
    rmSS_Total = sum((rand_Y - mean(rand_Y)).^2);
    randmodel_r_squared = 1 - (rmSS_Residual/ rmSS_Total);
    Random_modelR(j,1)=randmodel_r_squared;
    random_rvalue(j,1)=corr(rand_Y,random_y_pred);
    
end
save([data_dir,'ElasRandom_model_res.mat'],'Random_modelR','random_rvalue');

% calculate model significance
model_R2_pspin=(length(find(Random_modelR>real_model_Rsquare))+1)/(spins+1);
model_rvalue_pspin=(length(find(random_rvalue>real_rvalue))+1)/(spins+1);

save([data_dir,'Elasmodel_pspinres.mat'],'model_R2_pspin','model_rvalue_pspin');

%% plot observed and predicted MV brain maps
inor_Y=inormal(Y);
surf_plot(inor_Y,'bwocolorbar.mat')
inor_yped=inormal(yped);
surf_plot(inor_yped,'bwocolorbar.mat')

