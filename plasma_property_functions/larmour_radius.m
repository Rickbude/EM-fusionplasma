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
function [rho_L] = larmour_radius(T,m,q,B)
%LARMOUR_RADIUS Summary of this function goes here
%   Detailed explanation goes here
    v_th   = thermal_velocity(T,m);
    Omega_c  = cyclotron_frequency(q,B,m);
    rho_L  = v_th./Omega_c;
end

