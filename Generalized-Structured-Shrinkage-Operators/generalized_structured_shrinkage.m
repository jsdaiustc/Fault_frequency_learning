function [ x ] = generalized_structured_shrinkage( y, T, Neigh_Size, Method)
% This function calculate the generalized structured shrinkage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y:                          The input signal
% T:                          The K-sparsity threshold
% Neigh_Size:                 The size of the window
% Params: a struct contains all the parameters
%       Params.W_type:        The weight type: 'SESK', 'multi-scale PMI' or 'None'
%                             SESK: square envelope spectrum kurtosis (default)
%                             multi-scale PMI: multi-scale periodic modulation intensity
%                             None: the weight is inoperative
%       Params.Fs:            The sampling frequency of the simulation signal
%       Params.F_Interval:    The searching interval of the characteristic frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weight:                     The generated weight vector

% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2018.6
if nargin < 4
    Method = 'gausswin';
end
if nargin < 3
    Neigh_Size = 5;
end

y = y(:)';
L = length(y);
half = fix(Neigh_Size/2);
yy = [fliplr(y(2: 1 + half)), y, fliplr(y(end-half+1: end))];

Y = buffer(yy, Neigh_Size, Neigh_Size-1, 'nodelay');
Y = Y(:, 1 : L);  
Y = Y.^2;
switch Method.window
    case 'gausswin'
        W = window(@gausswin,Neigh_Size); 
    case 'triang'
        W = window(@triang,Neigh_Size); 
    case 'rectwin'
        W = window(@rectwin,Neigh_Size); 
    case 'hamming'
        W = window(@hamming,Neigh_Size); 
    otherwise
        warning('There is no such window.')       
end

W_norm = W / norm(W);
Y = bsxfun(@times, Y, W_norm);
Y_norm = sqrt(sum(Y, 1));
Y_norm = Y_norm(:)';
%% This is the shrinkage
% x = max(1 - (T./abs(Y_norm)).^2, 0) .* y;
Method.lambda = T*abs(y)./abs(Y_norm);
phi = Shrinkage(Method.SubName, Method);
x = phi(y(:));
x = x(:)';
end
