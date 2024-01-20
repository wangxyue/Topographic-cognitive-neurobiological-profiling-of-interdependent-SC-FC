%% Spintest
SphereSurf=cell(2,1);
SphereSurf{1}='Q1-Q6_RelatedParcellation210.L.sphere.32k_fs_LR.surf.gii';
SphereSurf{2}='Q1-Q6_RelatedParcellation210.R.sphere.32k_fs_LR.surf.gii';

LV=gifti('Glasser180_210P_L.func.gii'); % FS_LR 32k
RV=gifti('Glasser180_210P_R.func.gii');
NumHemi=size(LV.cdata, 1);


GlasserLR=[LV.cdata; RV.cdata];
LabelLR=zeros(size(GlasserLR));
LabelLR=GlasserLR;

Label=cell(2,1);
Label{1}=LabelLR(1:NumHemi, 1);
Label{2}=LabelLR(NumHemi+1:end, 1);

NumPerm=10000; 
for i=1:NumPerm
    OutFile=['PermLabel\Rand_' sprintf('%.5d', i)];
    fprintf('Performing %s\n', OutFile);
    GetRotateLabel(SphereSurf, Label, OutFile);
end
