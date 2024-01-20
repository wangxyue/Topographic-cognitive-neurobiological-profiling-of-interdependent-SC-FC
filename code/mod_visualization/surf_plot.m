function surf_plot(var,varargin)
%% node2txt
nVarargs = length(varargin);

gii1 = gifti('mod_visualization\plot_surface\Glasser180_210P_L.label.gii');
gii2 = gifti('mod_visualization\plot_surface\Glasser180_210P_R.label.gii');
ParcelLabel = double([gii1.cdata;gii2.cdata]);
Z= zeros(length(ParcelLabel),1);
for i = 1:360
    Z(ParcelLabel==i) = var(i,1); %把第i个脑区的指标值赋值给Z中对应位置处
end
save('mod_visualization\plot_surface\test.txt','Z','-ascii');
if nVarargs == 0
    BrainNet_MapCfg('mod_visualization\plot_surface\FSaverage_inflated_32K.nv','mod_visualization\plot_surface\test.txt','mod_visualization\plot_surface\Surf.jpg');    
else
    BrainNet_MapCfg('mod_visualization\plot_surface\FSaverage_inflated_32K.nv','mod_visualization\plot_surface\test.txt',strcat('mod_visualization\plot_surface\',varargin{1}));
end


