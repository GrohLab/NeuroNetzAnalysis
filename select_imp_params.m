function params_new = select_imp_params(params,I,S)
[imp,dnd]=sort(params(:,1),'ascend');
i = 1;
while (imp(i) < I && sum(params(:,1)) >= S) ||...
        params(dnd(i),3) < mean(params(:,3)) - std(params(:,3))
    params(dnd(i),:)=[];
    [imp,dnd] = sort(params(:,1),'ascend');
end
lowStd = params(:,3) < mean(params(:,3)) - std(params(:,3));
params_new = params(~lowStd,:);
params_new(:,1) = params_new(:,1)/sum(params_new(:,1));
end

