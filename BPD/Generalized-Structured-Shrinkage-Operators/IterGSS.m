function [ x ] = IterGSS(y, A, AH, normA, rho, K, Method, Nit, init)
% This function realizes iterative generalized structured shrinkage based
% on TQWT
%% Input %%%%%%%%%%
%   y      : noisy data
%   A, AH  : function handles for A and its transpose
%   normA  : the norm of the transformation
%   rho    : rho >= maximum eigenvalue of A'A
%   K      : K-sparsity parameters
%   Method : a structure with information of the method
%   Nit    : number of iterations
%   init   : the inital point of x
%% Output %%%%%%%%%%
%   x      : Recovery

% Reference: 'Sparsity-assisted Fault Feature Enhancement: Algorithm-aware versus Model-aware',
% IEEE Transactions on Instrumentation and Measurement, 2020
% https://zhaozhibin.github.io/
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2019.6

% initialization
AHy = AH(y);
if ~exist('init', 'var')
    init = AHy;
end
if ~exist('Nit', 'var')
    Nit = 50;
end
% cost = zeros(Nit , 1);
mu = 1 / rho;
x = init;
x_old = init;
Ax = A(x);
iter = 1;

while iter <= Nit
%     fprintf(['Iteration:' num2str(iter) '\n'])
    % forward operator
    % x = x - mu * (A^T(A(x) - y))
    tmp = x; 
    AHAx = AH(Ax); 
    for i = 1:numel(tmp) 
        tmp{i} = tmp{i} + mu * ( AHy{i} - AHAx{i} );  
    end
     
    % backward (thresholding)
    % x = soft(x, mu * lam);
    Temp = [];
    for i = 1:numel(x)
        Temp = [Temp ; tmp{i}(:) / normA(i)];
    end
    c = K_sparsity(Temp, K);
    for i = numel(x):-1:1
        if ~strcmp(Method.Name, 'WGL')
            Method.lambda = c * normA(i);
            Phi = Shrinkage(Method.Name, Method);
            x{i} = Phi(tmp{i});               
            
        else
            if i == numel(x)
                Size = Method.Initial_Size;
                x{i} = generalized_structured_shrinkage(tmp{i}, c * normA(i), Size ,  Method);
            else
                Length_Before = length(tmp{i+1});
                Length_Now = length(tmp{i});
                if Length_Before ~= Length_Now
                    Size = (Size - 1) * 2 + 1;
                end
                x{i} = generalized_structured_shrinkage(tmp{i}, c * normA(i), Size, Method);
            end
        end
    end    
 
    % reconstruction
    Ax = A(x);
    
    stop = 0;
    total = 0;
     for i = numel(x)
        stop = stop + sum((x{i} - x_old{i}).^2);
        total = total + sum((x{i}).^2);
    end   

    if sqrt(stop) / sqrt(total) < 1e-6
        break;
    end
    iter = iter + 1;
    x_old = x;
end

111;

