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
function [Sm] = Poynting_vector(E,B)
%POYNTING_FLUX Calculates the time-averaged "phasor" Poynting vector
%   Sm = 1/(2*mu0)*[ Em x Hm* ]
    mu0 = pi*4e-7;
    
    %E and H* phasor
    Em  = E;
    Hm  = conj(B)/mu0;
    
    Sm  = zeros(size(E));
    %poynting vector components (calculated by cross product Em x Hm* )
    Sm(:,1) = 1/2*(Em(:,2).*Hm(:,3) - Em(:,3).*Hm(:,2));    %Smx
    Sm(:,2) = 1/2*(Em(:,3).*Hm(:,1) - Em(:,1).*Hm(:,3));    %Smy
    Sm(:,3) = 1/2*(Em(:,1).*Hm(:,2) - Em(:,2).*Hm(:,1));    %Smz
end

