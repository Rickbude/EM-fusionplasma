function [Lambda] = Lambda_n(n,z)
%LAMBDA_N Summary of this function goes here
%   Detailed explanation goes here
    Lambda = besseli(n,z).*exp(-z);
end

