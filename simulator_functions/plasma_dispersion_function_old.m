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
%PLASMA_DISPERSION_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
    approach = 3;
    if(approach == 1)
        dz = 0.01;
        dt = 0.01;
        r  = 1e-14;
        z1 = -100:dz:zeta-r;
        t  = 0:dt:1;
        z2 = zeta+r*exp(1i*t*pi);    
        z3 = zeta+r:dz:100 ;



        out1 = exp(-z1.^2)./(z1-zeta);
        out2 = exp(-z2.^2);
        out3 = exp(-z3.^2)./(z3-zeta);
       % Z = sum(out1*dz) + -1i*pi*sum(out2*dt) +sum(out3*dz);
        Z = sum(out1*dz) + -1i*pi*exp(-zeta^2) +sum(out3*dz);
        Z = 1/sqrt(pi)*Z;
    elseif(approach == 2)
        if(abs(zeta)<6)
            dz = 0.0001;
            z  = 0:dz:abs(zeta);
            %S  = exp(-zeta^2)*sum(dz*exp(z.^2));
            S  = dawson(abs(zeta));
        else
            S = 1./(2*abs(zeta)) + 1./(4*abs(zeta).^3) + 3./(8*abs(zeta).^5);
        end
                
        if(zeta>0)
             Z  = 2*S  +1i*sqrt(pi)*exp(-zeta^2);
        else
             Z  = -2*S +1i*sqrt(pi)*exp(-zeta^2);
        end
        
    else
        %Z =  1i*sqrt(pi)*exp(-zeta.^2) - 2*dawson(zeta);
        %profiling show that the slowest factor here is the exp function
        %Speed gain can be made here
        
        S = dawson_approx(zeta);
        %S = dawson(zeta);
        
        Z =  1i*sqrt(pi)*exp(-zeta.^2) - 2*S;
    end
end

