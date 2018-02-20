function [params,threshold]=emforgmm(data,M,err,clin,varargin)
%% EMFORGMM Uses the expectation-maximitation algorithm fit M gaussians.
% emforgmm(DATA,M) returns the parameters for M gaussian models fitted for
% the DATA with random initialization.
% emforgmm(DATA,M,MU_1,SIG_1,...,MU_M,SIG_M) returns the parameters for M
% gaussian models fitted for the DATA with user defined initialization as
% MU_1, SIG_1,...,MU_M, SIG_M.
%% Function initialization
params = zeros(M,3);
data = data(:);
Maux = M;
MX = max(data)-std(data);
MN = min(data)+std(data);
data = sort(data(:),'ascend');
if isempty(varargin)
    params = random_init_params(M,MX,MN,std(data));
elseif length(varargin) == 2*M
    %% User defined initialization
    for l=1:M
        params(l,2) = varargin{l*2-1};
        params(l,3) = varargin{l*2};
    end
else
    error(['Initialization for the variables is not correct.',...
        ' Make sure that the mean and variance are set for all gaussians'])
end
% ALPHA
params(:,1) = 1/M * ones(M,1);
MAXREP = 513;
if exist('err','var')
    epsilon = err;
else
    epsilon = 1e-7;
end
N = length(data);
N_inv = 1/N;
sqerr = 1;
j=1;
gauss_l = zeros(M,N);
px_l = gauss_l;
pik = px_l;
L_old = 0;
h = waitbar(0,'Estimating GMM');
%% Starting EM for GMM
while j <= MAXREP && sqerr > epsilon
    params_old = params;
    for l = 1:M
        gauss_l(l,:) = evalgauss(data,params_old(l,2),params_old(l,3));
        px_l(l,:) = params_old(l,1)*gauss_l(l,:);
    end
    p_xhat = sum(px_l,1);
    if sum(isnan(p_xhat)) > 0
        break;
    end
    for l=1:M
        pik(l,:) = px_l(l,:)./p_xhat;
    end
    sum_pik = sum(pik,2);
    sum_pik_inv = 1./sum_pik;
    params(:,1) = N_inv*sum_pik;
    l=1;
    while l<=M
        params(l,2)=updatemu(data,pik(l,:),sum_pik_inv(l));
        params(l,3)=updatesig(data,params(l,2),pik(l,:),sum_pik_inv(l));
        if params(l,2) == 0  || params(l,3) == 0 ||...
                isnan(params(l,2)) || isnan(params(l,3)) ||...
                isinf(params(l,2)) || isinf(params(l,3))
%             If the estimation results in a zero, the GMM will be reduced
%             to M - 1.
            warning(['Deleting one component. Keeping ',num2str(M-1)])
            params(l,:)=[];
            params_old(l,:)=[];
            sum_pik_inv(l) = [];
            px_l(l,:) = [];
            pik(l,:) = [];
            M = M - 1;
            if M == 0
                params_old = random_init_params(Maux,MX,MN,std(data));
                params = params_old;
                disp('All parameters were deleted. Re-estimating...')
                M = Maux;
                l = 1;
                break;
            end
        else
            l = l + 1;
        end
    end
    L = sum(log(p_xhat));
    if isnan(L)
        break;
    end
    sqerr = abs(L_old-L)^2;
    L_old = L;
    if isnan(L_old )
        params_old = random_init_params(M,MX,MN,std(data));
        j = 1;
        disp('Calculation error! Restarting the process')
    end
    waitbar(sqerr\epsilon,h)
    j = j+1;
end
close(h)
likelihood = sum(p_xhat)*N_inv;
display(['Final likelihood: ',num2str(likelihood)])
% Selecting those parameters with a contibution bigger than 1% but
% constrainting that at least the 99% is kept.
if clin || likelihood < 0
    warning('Cleaning the parameter matrix.')
    I = 0.005;
    params = select_imp_params(params,I,1-I);
end
% plotGMMs(params,data)
if (~isempty(params) || sum(params(:,1)) ~= 0) && nargout == 2
    threshold = findthreshGMM(params, data);
    % threshold = 1;
else
    disp('Not a lucky day...')
    threshold = NaN;
end

end
% || params(l,2) == 0 
function params_init = random_init_params(M,MX,MN,DS)
params_init = zeros(M,3);
% WEIGHT
params_init(:,1) = 1/M;
% MEAN
params_init(:,2) = MN:(MX-MN)/M:MX-((MX-MN)/M);
% STANDARD DEVIATION
params_init(:,3) = rand(M,1)*DS;
end

function mu_new = updatemu(data,pik_mu,sum_pik_inv_mu)
% Example for a vector sum: 'taking advantage of the CPU architecture'
% suma1=0;
% suma2=0;
% suma3=0;
% suma4=0;
% for i=1:4:N-3
%     suma1 = suma1 + (data(i) * pik_mu(i));
%     suma2 = suma2 + (data(i+1) * pik_mu(i+1));
%     suma3 = suma3 + (data(i+2) * pik_mu(i+2));
%     suma4 = suma4 + (data(i+3) * pik_mu(i+3));
% end
% suma1 = suma1*sum_pik_inv_mu;
% suma2 = suma2*sum_pik_inv_mu;
% suma3 = suma3*sum_pik_inv_mu;
% suma4 = suma4*sum_pik_inv_mu;
% mu_new = suma1+suma2+suma3+suma4;
mu_new = (pik_mu*data)*sum_pik_inv_mu;
end

function sig_new = updatesig(data,mu_new,pik_sig,sum_pik_inv_sig)
sig_new = (pik_sig * ((data - mu_new).^2)) * sum_pik_inv_sig;
end