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
function [R,L,P,S,D] = cp_Stix_components(omega,B,m,n,q)
%CP_MODIFIED_STIX_COMPONENTS Summary of this function goes here
%   Detailed explanation goes here

    

    %plasma properties
    Omega   = cyclotron_frequency(q,B,m);  
    omega_p = plasma_frequency(n,q,m);
    
    %calculate the components
    
    nu    = imag(omega);
    omega = real(omega);
    
    
    R = 1-sum(omega_p.^2./(omega.*(omega+1i*nu+sign(q).*Omega)),2);
    L = 1-sum(omega_p.^2./(omega.*(omega+1i*nu-sign(q).*Omega)),2);
    P = 1-sum(omega_p.^2./(omega.*(omega+1i*nu               )),2);
    S = 1/2*(R+L);
    D = 1/2*(R-L);
end

