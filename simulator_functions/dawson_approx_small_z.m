%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019-2020 Eindhoven University of Technology.             %
%                                                                         %
% This code is free software, you can redistribute it and/or modify it    %
% under the terms of the GNU General Public License; either version 3.0   %
% of the License, or (at your option) any later version. See LICENSE.md   %
% for details, or see <https://www.gnu.org/licenses/>                     %
%                                                                         %
% This code is distributed in the hope that it will be useful, but        %
% WITHOUT ANY WARRANTY; without even the implied warranty of              %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    %
%                                                                         %
% Author:  Rick Budé       (r.h.s.bude@tue.nl)                            %
%                                                                         %
% Contact: Rick Budé       (r.h.s.bude@tue.nl)                            %
%          Jan van Dijk    (j.v.dijk@tue.nl)                              %
%          Roger Jaspers   (r.j.e.jaspers@tue.nl)                         %
%          Bart Smolders   (a.b.smolders@tue.nl)                          %
%          Dirk Van Eester (d.van.eester@fz-juelich.de)                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out] = dawson_approx_small_z(z,h,N)
%DAWSON_APPROX approximates the Dawson function for low(z<~30) z
%   % This approximation is discussed on:
    % http://www.ebyte.it/library/codesnippets/DawsonIntegralApproximations.html
    % Which in turn gets it from George B. Rybicki:
    % Dawson's Integral and the Sampling Theorem (1989)
    % Computers in Physics 3, page 85-87 
    % doi: https://doi.org/10.1063/1.4822832
    
    % The expansion runs from -N_odd < n < N_odd, for odd n
    % F(z) ~ 1/sqrt(pi) * sum_{n odd} 1/n*exp(-(z-nh)^2)
    % There is probably a - missing in eq (12) in Rybicki's paper
    
    %Depending on the value of z, more or less terms need to be included
    %e.g. for z<25: h = 0.4, N = 75 works very well.
    out = zeros(size(z));
       
    
    N_odd = N;
    %only consider odd N
    if(mod(N,2)==0)
       N_odd = N-1; 
    end
    
    %do the calculation. About 80% of this function is spent calculating
    %the exponents.
    % exp(-(z-nh)^2) = exp(-z^2)*exp(2*n*h*z)*exp(-n^2h^2)
    % As mentioned in Rybicki's paper:
    % - exp(-z^2)can be precomputed
    % - exp(-n^2*h^2) is cheap to compute 
    % - exp(2*h*n*z) is what takes ~80% of the time here
    % - exp(2*h*n*z) and exp(-2*h*n*z) are calculated simulataneously:
    % only calculate for positive n
    z_exp = exp(-z.^2);
    
    for n = 1:2:N_odd      
       z_part = exp(2*h*n*z);
       %"mirror" for negative n:
       %out = out + 1/n*exp(-(z-n*h).^2)
       out = out+1/n*(z_part - 1./z_part)*exp(-n^2*h^2);     
    end    
    out = 1/sqrt(pi)*z_exp.*out;
end

