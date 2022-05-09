function [ T ] = K_sparsity( x, k )
% This function find the threshold for K-Sparsity strategy

%% Input %%%%%%%%%%
%   x      : the input sequence
%   k      : the k-sparsity parameter
%% Output %%%%%%%%%%
%   T      : the generated shrinkage
% Author : Zhibin Zhao
% Place : Xi'an Jiaotong University
% Email : zhaozhibin@stu.xjtu.edu.cn
% Date : 2018.6
Descend = sort(abs(x), 'descend');
%% This is the percentage of the coefficients
% N = length(x);
% K = round(N * k);
K = k;
T = Descend(K);

end

