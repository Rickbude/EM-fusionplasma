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
    
    Bz = B(:,3);

    %load the Stix components
    [~,~,P,S,D]   = cp_Stix_components(omega,Bz,m,n,q);
    
    epsilon       = zeros(length(m),3,3);
  
    % magnetic field in z-direction
    epsilon(:,1,1) = S;
    epsilon(:,2,2) = S;
    epsilon(:,3,3) = P;
        
    epsilon(:,1,2) = -1i*D;
    epsilon(:,2,1) =  1i*D;        
    

end

