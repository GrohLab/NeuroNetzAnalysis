function plotGMMs(parameters,data)
data = sort(data(:),'ascend');
M=size(parameters,1);
for k=1:M
    %pik(k,:)=parameters(k,1)*evalgauss(data,parameters(k,2),...
    %    parameters(k,3));
    pik(k,:)=evalgauss(data,parameters(k,2),parameters(k,3));
    pik(k,:)=pik(k,:)*parameters(k,1);
end
p_xhat = sum(pik,1);
figure('name','Gaussian Mixture Model');
h=histogram(data,ceil(sqrt(length(data))),'Normalization','pdf','FaceColor',...
    [1,1,1]);hold on;
puk=plot(data,pik,'--','LineWidth',2.5);hold on;
pxh=plot(data,p_xhat,'r','LineWidth',3.0);
end