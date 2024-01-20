%% Neurocognitive flexibility association analysis
% define work dir
data_dir = 'G:\interdependent_SCFCnet\multi_modanalysis\analysis_data\';
load([data_dir,'pro_node.mat']);
load([data_dir,'multimod_detecres\r1w1para\cross_sub_MV.mat']);
load([data_dir,'SpinTest_result\finalres\Perm_360_MVResults.mat']);
node_flexibility=pro_node;
realR=corr(cross_sub_MV,node_flexibility,'type','pearson');
%% compute the significance
for perm=1:size(Perm_360_Results,2)
    perm_MV=Perm_360_Results(:,perm);
    r=corr(node_flexibility,perm_MV,'type','pearson');
    perm_R(1,perm)=r;
end
p_sig=(length(find(perm_R>realR))+1)/(perm+1);
%surf_plot(pro_node,'mycolor_bluered2.mat')
%% plot rvalue
h = histogram(perm_R, 80);
h.FaceColor = [0.3 0.3 0.3]; 
h.EdgeColor = [0.3 0.3 0.3];
set(gca, 'XTick', []);
set(gca, 'YTick', []);
box off
hold on
line([realR realR], [0 400], 'LineStyle', '--', 'Color', 'r');

%% four type nodes
low_nodeID=find(node_flexibility>=0 & node_flexibility<1);
moderate_nodeID=find(node_flexibility>=1 & node_flexibility<2);
good_nodeID=find(node_flexibility>=2 & node_flexibility<3);
high_nodeID=find(node_flexibility>=3);

MV1=cross_sub_MV(low_nodeID);
MV2=cross_sub_MV(moderate_nodeID);
MV3=cross_sub_MV(good_nodeID);
MV4=cross_sub_MV(high_nodeID);

g1=repmat({'MV1'},size(MV1));
g2=repmat({'MV2'},size(MV2));
g3=repmat({'MV3'},size(MV3));
g4=repmat({'MV4'},size(MV4));

data=[MV1;MV2;MV3;MV4];
group=[g1;g2;g3;g4];

% perform Kruskal-Wallis test
[p1, table1, stats1] = kruskalwallis(data,group);
% multiple comparison
[c,m,mh,gnames]=multcompare(stats1);
