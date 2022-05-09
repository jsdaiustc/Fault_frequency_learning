function [T] = Estimate_Period(Sig , Fn_N)
% Estimate the period of the signal
% Input:
%          Sig    : the input signal
%          Fn_N   : a vector which contains the period of each component (Fs / fc)
% Output:
%          T      : the estimated period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(Sig);
% hilbert 
% ht = abs(hilbert(Sig));
% difference hilbert
ht = abs(hilbert([0;diff(Sig)]));
ht_corr = xcorr(ht);
y = (ht_corr(N : end));
if norm(Fn_N) ~= 0
    Temp = -10 : 10;
    total_interval = [];
    for i = 1 : length(Fn_N)
        fc_interval = Temp + Fn_N(i);
        total_interval = [total_interval , fc_interval];
    end

    y_temp = zeros(N,1);
    y_temp(total_interval) = y(total_interval);
    [~ , T] = max(y_temp);
    T=T-1; % NOTE modified !!!  
    
else
% the below is directly search the period, but it is not suitable to high
% noise
    k = 1;
    % delete the value near the r(0)
    for i = 1 : N-1
        if y(i+1) <= y(i)
            y(i) = 0;
        else
            break;
        end
    end
    % search the maximum index
    y_diff = diff(y);
    for i = 1 : length(y_diff)-1
        if y_diff(i)>=0 && y_diff(i+1)<= 0
            index(k) = i;
            k = k + 1;
        end
    end
    y_final = zeros(N,1);
    y_final(index) = y(index);
    [~ , T] = max(y_final);
%     [~ , T] = max(y);
end
end

