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

function [backward,central,forward] = fd_coefficients(N)
%FD_COEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here
    
    n         = 2; %accuracy. This is currently still fixed at 2
    Nmax      = 2*floor((N+1)/2)-1+n;
    central   = zeros(N+1,Nmax);
    backward  = zeros(N+1,N+n);
    forward   = zeros(N+1,N+n);
    for m = 1:N
        
        %%%%%%%%%%%%%%%% FORWARD DIFFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        ncoeff = m+n;
        s      = 0:ncoeff-1;
        coeff = fd_coefficients_arbitrary_stencil(m,s);
        forward(m+1,s+1) = coeff;            
                
        %%%%%%%%%%%%%%%% BACKWARD DIFFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ncoeff= m+n;
        s     = -(ncoeff-1):0;
        coeff = fd_coefficients_arbitrary_stencil(m,s);
        backward(m+1,N+n+s) = coeff;               

        %%%%%%%%%%%%%%%% CENTRAL DIFFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %aim for nth order accuracy
        %ncoeff = 2*floor((m+1)/2)-1+n;
        %aim for maximum accuracy with N stencil points
        ncoeff = Nmax;
        p      = floor(ncoeff/2);        
        s      = -p:p;        
        coeff = fd_coefficients_arbitrary_stencil(m,s);        
        central(m+1,ceil(Nmax/2)+[-p:p]) = coeff;        
        
    end
    
    %also set numbers for 0-th order derivative
    central(1,ceil(Nmax/2)) = 1;
    backward(1,N+n) = 1;
    forward (1,  1) = 1;
    
    
    
end

