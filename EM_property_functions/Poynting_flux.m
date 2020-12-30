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
function [Flux] = Poynting_flux(r,E,B,curve_k_parr)
%POYNTING_FLUX Poynting flux, currently only for (-?) x-direction
%   Detailed explanation goes here
    Sm = Poynting_vector(E,B);
    
    Flux      = zeros(size(Sm));
    if(curve_k_parr)
        Flux(:,1) = -Sm(:,1)*2*pi.*r;
    else
        Flux(:,1) = -Sm(:,1);
    end
    
    %TODO: x and y component are not calculated yet
end

