function [ phi ] = penalty( pen, lambda )
% This function provides different penalties.
%% Input %%%%%%%%%%
%   pen      : the name of penalty
%   lambda   : the trade-off parameter (can be a vector)
%% Output %%%%%%%%%%
%   phi      : the function handle of the penalties
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2019.6
lambda = lambda(:);
switch pen
    case 'L1'
        phi = @(x) lambda .* abs(x);
    case 'L0'
        phi = @(x) lambda .* (abs(x) > 0);
    case 'Lp'
        phi = @(x , p) lambda .* abs(x) .^ p;
    case 'SCAD'
        phi = @(x , a) (lambda .* abs(x)) .* (abs(x) < lambda) ...
            + (2*a*lambda.*abs(x)-x.^2-lambda.^2) / (2*(a-1)) .* (abs(x) >= lambda & abs(x) < a*lambda) ...
            + ((a+1) * lambda.^2 / 2) .* (abs(x) >= a*lambda);
    case 'MC'        
        phi = @(x , gamma) lambda .* (abs(x) - x.^2  ./ (2*gamma*lambda)) .* (abs(x) <= gamma*lambda) ...
            + gamma*lambda.^2/2 .* (abs(x) > gamma*lambda);
    otherwise
        error(' penalty must be ''L1'', ''L0'', ''Lp'', ''SCAD'', ''MC'' ')
end

end

