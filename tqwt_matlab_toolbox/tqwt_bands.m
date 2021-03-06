function y = tqwt_bands(b,w,Q,r,J,N)
% y = tqwt_bands(b,w,Q,r,N)
% Reconstruction from subsets of TQWT subbands
% INPUT
%   b - cell array (length K) of subbands
%   w - TQWT coefficients
%   Q, r, J - TQWT parameters
%   N - length of signal.
% OUTPUT
%   y - array of reconstructed signals (one per row)
%
% % Example



K = length(b);

y = zeros(K,N);

wz = cell(1,J+1);
for j = 1:J+1
    wz{j} = zeros(size(w{j}));
end

for k = 1:K
    bk = b{k};
    for i = 1:length(bk)
        m = bk(i);
        w2 = wz;
        w2{m} = w{m};
        y(k,:) = y(k,:) + itqwt_radix2(w2,Q,r,N);
    end
end



