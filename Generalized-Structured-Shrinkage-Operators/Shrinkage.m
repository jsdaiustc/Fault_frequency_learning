function [ shrink ] = Shrinkage( pen, params )
% This function provides the shrinkage operators of different penalties.
%% Input %%%%%%%%%%
%   pen      : the name of penalty
%   params   : the structure of the parameters
%         params.lambda  : the trade-off parameter (can be a vector)
%         params.p       : the parameter in lp
%         params.SCAD    : the parameter in SCAD, should be larger than 2
%         params.MC      : the parameter in MC, should be larger than 1
%% Output %%%%%%%%%%
%   shrink      : the function handle of the shrinkage operators
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2019.6

lambda = params.lambda(:);
switch pen
    case 'L1'
        shrink = @(x) sign(x) .* max(abs(x)-lambda, 0);
    case 'L0'
        shrink = @(x) 0 .* (abs(x) < sqrt(2*lambda)) ...
            + x .* (abs(x) >= sqrt(2*lambda));
    case 'Lp' % here we use p = 1/2 and 2/3 
        shrink = @(x) P_value(x , params.p, lambda);
    case 'SCAD'
        shrink = @(x) sign(x).*max((abs(x)-lambda),0) .* (abs(x) <= 2*lambda) ...
            + ((params.a-1)*x-sign(x)*params.a*lambda)/(params.a-2) .* (abs(x) > 2*lambda & abs(x) <= params.a*lambda) ...
            + x .* (abs(x) > params.a*lambda);
    case 'MC'   % mu = lambda * gamma  for firm  
        shrink = @(x) 0 .* (abs(x) <= lambda) ...
            + sign(x).*(abs(x)-lambda)/(1-1/params.gamma) .* (abs(x) > lambda & abs(x) <= params.gamma*lambda)...
            + x .* (abs(x) > params.gamma*lambda);
    otherwise
        error(' Shrinkage must be ''L1'', ''L0'', ''Lp'', ''SCAD'', ''MC'' ')
end
end

function x = P_value(x, p, lambda)
tmp = x(:);
n = length(x);
switch p
    case 1/2
        Index = find(abs(tmp) > (54)^(1/3)/4*lambda^(2/3));
        x = zeros(n,1);
        x(Index) = real(2/3*tmp(Index).*(1 + cos(2/3*(pi - acos(lambda/8.*((abs(tmp(Index))/3).^(-1.5)))))));
    case 2/3
        Index = find(abs(tmp)>(2/3)*(3*(lambda)^3)^(1/4));
        x = zeros(n,1);
        x(Index) = 2/sqrt(3)*(lambda).^(0.25).*(cosh(acosh(27/16*tmp(Index).^2*(lambda).^(-3/2))./3)).^(0.5);
        x(Index) = sign(tmp(Index)).*real(((abs(x(Index))+sqrt(2*abs(tmp(Index))./abs(x(Index))-abs(x(Index)).^2))/2).^3);
end
end