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

function [epsilon] = cp_dielectric_tensor(omega,B,m,n,q)
%CP_DIELECTRIC_TENSOR returns the cold plasma dielectric tensor
%   Detailed explanation goes here
    
    Bx = B(:,1);
    By = B(:,2);
    Bz = B(:,3);
    
    %plasma properties
    %gyrotron frequencies
    Omega   = sign(q).*cyclotron_frequency(q,B ,m); %x+y+z
    Omega_x = sign(q).*cyclotron_frequency(q,Bx,m); %x-only
    Omega_y = sign(q).*cyclotron_frequency(q,By,m); %y-only
    Omega_z = sign(q).*cyclotron_frequency(q,Bz,m); %z-only
        
    %plasma frequency
    omega_p = plasma_frequency(n,q,m);
    
    %construct the 
    epsilon      = zeros(length(B),3,3);
    commonfact   = omega_p.^2./(omega.^2.*(Omega.^2-omega.^2));
    
    epsilon(:,1,1) = 1-sum(commonfact.*(Omega_x.*Omega_x + 1i*omega.*1i*omega),2);
    epsilon(:,1,2) =  -sum(commonfact.*(Omega_x.*Omega_y + 1i*omega.*-Omega_z),2);
    epsilon(:,1,3) =  -sum(commonfact.*(Omega_x.*Omega_z + 1i*omega.* Omega_y),2);
    
    
    epsilon(:,2,1) =  -sum(commonfact.*(Omega_y.*Omega_x + 1i*omega.* Omega_z),2);
    epsilon(:,2,2) = 1-sum(commonfact.*(Omega_y.*Omega_y + 1i*omega.*1i*omega),2);
    epsilon(:,2,3) =  -sum(commonfact.*(Omega_y.*Omega_z + 1i*omega.*-Omega_x),2);
    
    
    epsilon(:,3,1) =  -sum(commonfact.*(Omega_z.*Omega_x + 1i*omega.*-Omega_y),2);
    epsilon(:,3,2) =  -sum(commonfact.*(Omega_z.*Omega_y + 1i*omega.*Omega_x ),2);
    epsilon(:,3,3) = 1-sum(commonfact.*(Omega_z.*Omega_z + 1i*omega.*1i*omega),2);
    
end

