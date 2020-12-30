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

function [epsilon] = hp_dielectric_tensor_swanson(N_orders,omega,k,B,T,q,n,m)
%hp_dielectric_tensor_swanson calculates the hot plasma dielectric tensor
%as defined in Swanson

    N         = max(size(k,1),size(n,1));
    N_species = size(m,2);

    %The result is currently returned as an Nx(3x3) array. So a dielectric
    %tensor for each gridpoint, or for each wavenumber. Perhaps Nx9 is a
    %bit more readable, and more efficient
    epsilon = zeros(N,3,3);
    %pass the arguments, but load the tensor per species and then sum the
    %contributions
    for species = 1:N_species
        %extract temperature, charge, density, mass per species
        Tj = T(:,species);
        qj = q(:,species);
        nj = n(:,species);
        mj = m(:,species);
        %calculate the contribution from this single species
        epsilon     = epsilon + hp_dielectric_tensor_swanson_per_species(N_orders,omega,k,B,Tj,qj,nj,mj);
    end    
    %add the unit dyad to the dielectric tensor
    epsilon(:,1,1) = epsilon(:,1,1) + 1;
    epsilon(:,2,2) = epsilon(:,2,2) + 1;
    epsilon(:,3,3) = epsilon(:,3,3) + 1;    
end

