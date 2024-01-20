warning off

for i_sub = 1:60
temp_sub = MTC_Power264_Adult{i_sub,1};
%get timepoints window size=100s ,step=1TR
for k= 0:180  
     for j=(1+k):(50+k)
         temp1(j-k,:) = temp_sub(j,:); 
     end
     temp(k+1,1) = mat2cell(temp1);
end
FCPower1TR= gretna_fc_pearson(temp);

%稀疏矩阵
for i=1:181
FC_pos = FCPower1TR{i,1}.*(FCPower1TR{i,1}>0);      %只取正边
FC_bin = gretna_R2b(FC_pos,'s',0.1);    %将矩阵稀疏化，稀疏到0.1，此为二值矩阵 
FC_wei = FC_bin.*FC_pos; 
subwei1TR{i,1}=FC_wei;
end

gamma = 1; omega = 1; 
N=length(subwei1TR{1});
T=length(subwei1TR);
[B,mm] = multiord(subwei1TR,gamma,omega); 
[S,Q,n_it] = iterated_genlouvain(B,10000,0,0);
Q=Q/mm;
S = reshape(S, N, T);
% [S_sorted,s]=sort_ordinal(S);
% figure;
% imagesc(S_sorted);
modular_result(i_sub,1) = mat2cell(S);
all_Q(i_sub,1)=Q;
display(['subject ',num2str(i_sub),' is OK!']);
end

% %统计个数
% most_number=mode(S,2);

% %统计个数以及比例
% most_number=mode(S,2);
% b=[];
% [m n]=size(S);
% for i = 1 : m
% [k l]=mode(S(i,:));  %k为出现最多的数，l为出现最多的次数
% %b=[b;k l];
% b=[b;l];
% end
% b=b/181;
% 
% %看连通性
% z=[];
% for i=1:150
%   [m,n]=components(sparse(subwei08{i,1}));
%   z=[z;max(n)];
% end

% % single layer modular
% layer = 120;
% FC_gp = FCPower_60s{layer,1};
% Ci = S(:,layer);
% [~,Ind]= sort(Ci,'ascend');          %将分得的模块序号排序
% figure
% imagesc(FC_gp(Ind,Ind));
% unique(S(:,layer));
