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
function [Jant] = antenna_current(I_ant,r,r_a)
%ANTENNA_CURRENT give antenna current. Antenna sheath is assumed
%such that J_ant = delta(r-r_a)*(Jy y\hat + Jz z\hat)

    N        = length(r);
    h        = abs(max(r)-min(r))/(N-1);
    
    Jantx    = zeros(N,1);
    Janty    = zeros(N,1);
    Jantz    = zeros(N,1);
    
    [~,ant_ind]    = min(abs(r-r_a));
    
    I_ant_y  = I_ant(2);
    I_ant_z  = I_ant(3);
    
    Janty(ant_ind) = I_ant_y/h;
    Jantz(ant_ind) = I_ant_z/h;

    Jant     = [Jantx Janty Jantz];
    


end

