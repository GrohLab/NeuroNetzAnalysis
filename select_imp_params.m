function params_new = select_imp_params(params,I,S)
[imp,dnd]=sort(params(:,1),'ascend');
i = 1;
while (imp(i) < I && sum(params(:,1)) >= S) || params(dnd(i),3) < 1e-3
    params(dnd(i),:)=[];
    [imp,dnd] = sort(params(:,1),'ascend');
end
params_new = params;
params_new(:,1) = params_new(:,1)/sum(params_new(:,1));
end

