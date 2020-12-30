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
N               = 2001;         %number of grid points in k-space
N_high_res      = 50001;        %for generating pretty pictures
N_harmmax       = 4;            %number of terms to include in the dielectric tensor sum
iunityflux      = false;        %normalize flux, such that the maximum flux (at antenna) = 1
Nk              = 101;          %number of k-space sample points. N/20 is a reasonable choice in most cases
Np              = 20;           %polynomial order
zeta            = 2.2;         %k_perp*rho_L number for throughout domain
weak_form       = true;         %collocation model: use weak form formulation
lagrange_mult   = true;         %set boundary conditions as Lagrange multiplier type
cold_plasma     = false;        %use cold plasma tensor

%machine dimensions
a               = 2.20;         %device minor radius
ap              = 1.80;         %plasma minor radius
r0              = 6.20;         %major radius
x_a             = 1.90;         %antenna position relative to machine axis
x_wall          = 1.90;         %radial position of wall

%excitation
I_ant_y         = 1;            %antenna current (y)
I_ant_z         = 0;            %antenna current (z)
f               = 54e6;       %real excitation frequency
nuom            = 0e-3;         %collision frequency relative to omega
ntor            = 33 ;          %toroidal mode number (kz = ntor/r)
ky              = 0;            %perpendicular k component
curve_k_parr    = true;        %make k_parr dependent on r as n_tor/r
aperature       = true;         %use aperature-type excitation (1) or current sheet (0)


%plasma parameters
B0              = 5.30;         %magnetic field on axis
n0              = 6e19;         %core density[m^-3]
n0_edge         = 1e19;         %edge density, if a parabolic profile is used
alpha_n         = 0.15;         %density profile parameter if a parabolic profile is used
lambda_n        = 0.05;         %tune exponential decay in edge (smaller = faster decay)
T0              = 17e3;         %core temp (eV)
T0_edge         = 100;          %edge temperature (eV)
alpha_T         = 2.0;          %profile shpae parameter for temperature if parabolic profiles are used
lambda_T        = 0.04;         %tune exponential decay in edge (smaller = faster decay)
iatanprofile    = false;        %use "advanced density profile"
fracap          = 0.01;         %tune "advanced density profile"
                    
%ions in plasma:  He3   D       T          He4t       He4h    Be
n_species       = 6;            %number of ion species
A_ion           = [3       2       3       4        4            9];  %mass number of ions
Z_ion           = [2       1       1       2        2            4];  %charge of ions
conc_ion        = [0.04    0.44    0.44    0.06     0.01      0.01];  %concentration relative to ne
T_ion           = [T0      T0      T0      T0       600e3       T0];  %ion temperature.
alpha_T_ion     = [alpha_T alpha_T alpha_T alpha_T  0.1    alpha_T];
alpha_n_ion     = [alpha_n alpha_n alpha_n alpha_T  1.0    alpha_n];



