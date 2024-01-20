warning off
%to produce subject's modular across all windows
for i_sub = 1:506
    temp=modular_result{i_sub,1};
    for i=1:264
        for j=1:264
            a=temp(i,:)-temp(j,:);
            b=sum(a(:)==0);
            A(i,j)=b;
        end   
    end
    A = A-diag(diag(A));
    A = A/181;
    A(A<0.5) = 0; 
    temp1=mat2cell(A);
    
    gamma = 1; omega = 1; 
    N=length(temp1{1});
    T=length(temp1);
    B = multiord(temp1,gamma,omega); 
    [S,Q,n_it] = iterated_genlouvain(B,10000,0,0);
    S = reshape(S, N, T);
%     [S_sorted,s]=sort_ordinal(S);
%     figure;
%     imagesc(S_sorted);
    single_sub_modular_result(:,i_sub) = S;
    display(['subject ',num2str(i_sub),' is OK!']);
end


