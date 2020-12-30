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

%simulation parameters
N               = 5001;         %number of grid points in k-space
iunityflux      = false;        %normalize flux, such that the maximum flux (at antenna) = 1

%machine dimensions
a               = 2.10;         %device minor radius
ap              = 1.80;         %plasma minor radius
r0              = 6.20;         %major radius
x_a             = 1.98;         %antenna position relative to machine axis

%excitation
I_ant_y         = 1;            %antenna current (y)
I_ant_z         = 0;            %antenna current (z)
f               = 53e6;         %real excitation frequency
nuom            = 0e-3;         %collision frequency relative to omega
ntor            = 26 ;          %toroidal mode number (kz = ntor/r)
ky              = 0;            %perpendicular k component
curve_k_parr    = true;        %make k_parr dependent on r as n_tor/r

%plasma parameters
B0              = 5.30;         %magnetic field on axis
n0              = 1e20;         %core density[m^-3]
n0_edge         = 1e19;         %edge density, if a parabolic profile is used
alpha_n         = 1;            %density profile parameter if a parabolic profile is used
lambda_n        = 0.05;          %tune exponential decay in edge (smaller = faster decay)
T0              = 9e3;          %core temp (eV)
T0_edge         = 100;          %edge temperature (eV)
alpha_T         = 1.5;          %profile shpae parameter for temperature if parabolic profiles are used
lambda_T        = 0.20;         %tune exponential decay in edge (smaller = faster decay)
iatanprofile    = false;        %use "advanced density profile"
fracap          = 0.01;         %tune "advanced density profile"
                    
%ions in plasma:   H    D
n_species       = 3;            %number of ion species
A_ion           = [2    3     3   ];  %mass number of ions
Z_ion           = [1    1     2   ];  %charge of ions
conc_ion        = [0.49 0.49  0.02];  %concentration relative to ne






