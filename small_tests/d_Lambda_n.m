function [d_Lambda] = d_Lambda_n(K,n,z)
%D_LAMBDA_N Summary of this function goes here
%   K: order of derivative
%   n: "base" order of bessel functino
%   z: bessel function argument
    
    d_Lambda = zeros(size(z));
    for k = -K:K
        factor = (-1)^(k+K)*nchoosek(2*K,k+K);
        d_Lambda = d_Lambda + factor.*Lambda_n(n+k,z);
    end
    d_Lambda = d_Lambda / (2^K);
end

