%画出Q随run的次数增多的变化图
Q_max=[];
for run= 1:180
    for times =1:10
       temp = randsample(200,run,false);
       Q_max(run,times)=max(time_Q(temp,:));
       temp=[];
    end
    temp = [];
end
Q_max=Q_max';
Average=mean(Q_max);
Variance=var(Q_max);
%errorbar(1:180,Average(:,1:180),Variance(:,1:180))    %函数调用格式 errorbar(A,B,X)
plot(1:10:180,Average(:,1:10:180))
xlabel('run times');ylabel('maxQ');

%看多次平均的MV map的相似度
for i=1:200
    modular_result_transport_time{i,1}=modular_result_time{i,1}';
end
for run = 1:200
   temp = modular_result_transport_time{run,1};
   MV_temp = scaled_inclusivity_wei(temp,3);
   MV(run,1)=sum(MV_temp)/264;
   MV_all(run,:) = MV_temp';
   display(['run ',num2str(run),' is OK!']);
end

k=1;
MV_mean=[];
for run= 2:10:180    
    temp = randsample(200,run,false);
    MV_mean(k,:) = mean(MV_all(temp,:));
    temp=[];  
    k=k+1;
end

for i=1:17
    MV_Map_run_similar(i,1)= corr(MV_mean(i,:)',MV_mean(i+1,:)','type','Pearson');
end
plot(12:10:180,MV_Map_run_similar);
xlabel('run times');
ylabel('MV map similar');
