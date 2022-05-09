function [x , cost] = AdaESPGL( Sig, Params)
% This function realizes adaptive enhanced sparse period-group
% lasso (AdaESPGL) algorithm

% Input :
%       Sig : the input signal (1D)
% Params: a struct contains all the parameters
%       Params.N1:            The samples of one impulse
%       Params.M:             The number of periods
%       Params.Fn_N:          A vector which contains the period of each component (Fs / fc)
%       Params.lam:           The parameter related to sparsity across groups
%       Params.mu:            The parameter related to sparsity within groups
%       Params.pen:           The penalty function
%       Params.rho:           The degree of nonconvex
%       Params.Nit:           Number of iterations (default: 100)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       x:                    The denoised signal
%       cost:                 The history of the cost function
% 
% Output : 
%       x : the denosied signal
%       cost : the cost function 
%
% Reference: 'Enhanced sparse period-group lasso for bearing fault
% diagnosis', IEEE Transactions on Industrial Electronics, 2018
% HomePages: https://zhaozhibin.github.io/
% Author   : Zhibin Zhao
% Place    : Xi'an Jiaotong University
% Email    : zhibinzhao1993@gmail.com
% Date     : 2017.10
N1   = Params.N1;
M    = Params.M;
Fn_N = Params.Fn_N;
lam  = Params.lam;
mu   = Params.mu;
pen  = Params.pen;
rho  = Params.rho;
Nit  = Params.Nit;


% Fn_N=0;

% N1=4;
% M=6;

Sig = Sig(:);
x = Sig; % this initialization will avoid the zero at the denominator
cost = zeros(1 , Nit);

N_orig=length(x);

for iter = 1 : Nit
    P = Estimate_Period(x , Fn_N);
    P=min(round(length(x)/10),P);
 
    N0 = P - N1;
    b = binaryblock(1 , N0 , N1 , M );  % of size K*1  or   (N0+N1)*M   see EQ (21)
    N = sum(b);     
    a_max = 1/(lam*N);
    a_coefficient = rho * a_max;
    [phi , wfun] = pen_fun(a_coefficient , pen);
    norm_bx = sqrt( conv(abs(x).^2, b, 'full') );
    w = 1 ./ (abs(x) + eps);
    cost(iter) = 0.5 * sum(abs(x-Sig).^2) + lam * (sum(phi(norm_bx)) + mu * N * sum(abs(x) .* w));
    phi_derive_d_u = 1 ./ (wfun(norm_bx) + eps); % P2147 , r式中b_j右边的值， 注意wfun的定义是 按倒数定义的
    r = conv(phi_derive_d_u , b , 'valid');
    x = Sig ./ (1 + lam * r + eps);
    T = lam * mu * N ./ (1 + lam * r + eps) .* w;
    x_old=x;
    x = softth(x , T);    
    
   %% add

%     percent=0.95;
%     ind_zeros=find(x==0);
%     if length(ind_zeros)>N_orig*percent
%        x=zeros(N_orig,1);
%        [~, sort_ind]=sort(abs(x_old),'descend') ;
%        x(sort_ind(1:round(N_orig*(1-percent))))=x_old(sort_ind(1:round(N_orig*(1-percent))));   
%        break;
%     end
        
        
    
    
    
end

end




