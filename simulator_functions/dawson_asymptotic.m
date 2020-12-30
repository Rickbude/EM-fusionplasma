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
function [S] = dawson_asymptotic(zeta,N)
%DAWSON_ASYMPTOTIC dawson function for real zeta
%   Asymptotic expansion for the Dawson function
%   See Swanson, Plasma Waves (2nd ed.), Appendix B, eq. B.5
%   Use N to increase or decrease the number of terms included. 

%   For real zeta, the following values are appropriate
%   N=2: |zeta| > 1000
%   N=4: |zeta| > 100
%   N=6: |zeta| > 30
%   N=8: |zeta| > 15
    
    %Calculate inverse of zeta 
    zeta_inv        = 1./zeta;  
    n = 1:N;   
    
    zeta_inv_mat    = zeta_inv.^(2*n-1);
        
    coeff_top       = fact2(2*n-3);
    coeff_bottom    = 2.^n;
    coeff           = coeff_top./coeff_bottom;   
    
    S               = sum(coeff.*zeta_inv_mat,2);    
end

