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
function [B] = magnetic_field_1_5D(r,k,omega,E)
%CP_CALC_MAGNETIC_FIELD_1_5D Calculate the magnetic field using the
% Maxwell-Faraday equation.
%   Plane wave expansions in y and z are assumed 
%   Derivative with respect to x is determined using a central finite
%   difference 

    %electric field components
    Ex = E(:,1);
    Ey = E(:,2);
    Ez = E(:,3);
    
    %used k-components
    ky = k(:,2);
    kz = k(:,3);
    
    N  = length(E);
    h  = (max(r) - min(r))/(N-1);
    
    %derivative of E with respect to x (1st order central difference)
    %TODO: only fixed grid spacing is allowed
    %TODO: move this to a seperate function if it is used more often
    dE            = zeros(size(E));
    dE(2:N-1,:)   = (E(3:N,:)-E(1:N-2,:))/(2*h);
    %forward difference at leftmost node
    dE(1,:)       = (E(2,:)-E(1,:)  )/h;
    %backward difference at rightmost node
    dE(N,:)       = (E(N,:)-E(N-1,:))/h;
    
    dEy           = dE(:,2);
    dEz           = dE(:,3);    
    
    %calculate the curl. If k_parr_curve is active, kz(r) = kz/r
    curl_x = 1i*ky.*Ez - 1i*kz.*Ey;
    curl_y = 1i*kz.*Ex - dEz -1./r.*Ez;
    curl_z = dEy       - 1i*ky.*Ex;
    
    B = 1/(1i*real(omega))*[curl_x curl_y curl_z];
    
end

