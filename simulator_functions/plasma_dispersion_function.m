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

function [Z] = plasma_dispersion_function(zeta)
%PLASMA_DISPERSION_FUNCTION Calculate the plasma dispersion function
% The plasma dispersion function Z is related to the Dawson integral S via
% Z = i*sqrt(pi) * exp(-zeta^2) - 2S
% Or, it can be related to the error function 
% ( see Fried & Conte 1961: the plasma dispersion function ) 
% Z = i*sqrt(pi) * exp(-zeta^2) * (1 +  erf(i*z))
%   = i*sqrt(pi) * exp(-zeta^2) * (1 + i*erfi(z))

% The dawson integral is approximated, because it is a quite slow function
% to calculate. The approximations are not general, and are meant for cases
% with small imag(zeta)
            
    %S = dawson(zeta);
    S = dawson_approx_old(zeta);
    Z = 1i.*sqrt(pi).*exp(-zeta.^2) - 2*S;    
    
end

