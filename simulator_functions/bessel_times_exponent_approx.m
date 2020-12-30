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
function [approx] = bessel_times_exponent_approx(a,z)
%BESSEL_TIMES_EXPONENT_APPROX Approximates exp(-z)*besseli(a,z)
%   Two approximations, one for z < z_th, one for z > z_th
%   Detailed explanation goes here
    

    approx_low_z = 0;
    if(a<50)
        z_threshold = 100;
    else
        z_threshold = 500;
    end
    z_low    = z(abs(z)<z_threshold );
    z_high   = z(abs(z)>=z_threshold);

    %low z (z<100)
    
    if(approx_low_z)  
        approx_low = zeros(size(z_low));
        N = max(5*a,50);
        factors = -cumsum(log(1:N+a+1));
        for m = 0:N        
            if(m==0)
                if(a==0)
                    factor  = 0;
                else
                    factor  = factors(m+a);    
                end
            else
                factor  = factors(m+a) + factors(m);
            end
            approx_low = approx_low + exp(factor + (2*m+a)*log(z_low/2)-z_low);
        end
    else
       %Do not approximate the bessel function for low z
       approx_low = besseli(a,z_low).*exp(-z_low); 
    end
    
    %Approximate the exp(-z)*besseli(a,z) for high z
    %The Bessel function I(a,z) for fixed a and large |z| is given in 
    %Abramowitz and Stegun page 377, eq. 9.7.1, (|arg z|<pi/2)
    %http://people.math.sfu.ca/~cbm/aands/page_377.htm
    %I(a,z) ~ exp(z)/sqrt(2*pi*z)*{
    %                              + 1 
    %              - (mu-1)/(1!(8z)^1)  
    %        + (mu-9)(mu-1)/(2!(2z)^2)
    % - (mu-25)(mu-9)(mu-1)/(3!(8z)^3)
    % + ...
    %} = exp(z)/sqrt(2*pi*z)*{a0/x^0 + a1/x^1 + a2/x^2 + ...}  
    %With a0 = 1 and a_k = -a_{k-1} * (mu - (2n-1)^2)/(8n)
    %And with mu = 4a^2
    %The leading exp(z) in this approximation of course cancels out with
    %the exp(-z) from the required "composite function" exp(-z)*I(a,z)
        
    %Take 50 terms into account
    N  = 50;
    mu = 4*a^2;
    approx_high = ones(size(z_high));
    xn     = ones(size(z_high));
    ak     =1; 
    
    for n = 1:N  
        %update the coefficient from the previous generation
        ak     = -ak*(mu - (2*n-1)^2)/(n*8);
        xn     = xn./z_high;
        if(sum(isnan(ak*xn))==0)
            approx_high = approx_high+ak*xn;
        else
            break;
        end
    end

    approx_high = 1./(sqrt(2*pi*z_high)).*approx_high;

    
    approx                  = zeros(size(z));
    approx(abs(z)<z_threshold)   = approx_low;
    approx(abs(z)>=z_threshold)  = approx_high;
   

end

