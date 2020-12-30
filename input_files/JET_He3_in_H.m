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
N               = ceil(2000/18)*18+1; %number of grid points in k-space
N_high_res      = 10001;        %for generating pretty pictures
N_harmmax       = 4;            %number of terms to include in the dielectric tensor sum
iunityflux      = false;        %normalize flux, such that the maximum flux (at antenna) = 1
Nk              = 200;          %number of k-space sample points. N/20 is a reasonable choice in most cases
Np              = 8;            %polynomial order
zeta            = 1.2;          %k_perp*rho_L number for throughout domain
weak_form       = true;         %collocation model: use weak form formulation
lagrange_mult   = true;         %set boundary conditions as Lagrange multiplier type
cold_plasma     = false;        %use cold plasma tensor

%machine dimensions
as              = 1.35;         %simulation domain
a = as;
ap              = 0.95;         %plasma minor radius
r0              = 2.97;         %major radius
x_ant           = 1.05;         %antenna position relative to machine axis
x_wall          = 1.05;         %radial position of wall

%excitation
I_ant_y         = 1;            %antenna current (y)
I_ant_z         = 0;            %antenna current (z)
f               = 34.6e6;         %real excitation frequency
nuom            = 1e-3;         %collision frequency relative to omega
ntor            = 27 ;          %toroidal mode number (kz = ntor/r)
ky              = 0;            %perpendicular k component
curve_k_parr    = true;         %make k_parr dependent on r as ntor/r
aperature       = true;         %use aperature-type excitation (1) or current sheet (0)

%plasma parameters
B0              = 3.45;         %magnetic field on axis
n0              = 7e19;         %core density[m^-3]
n0_edge         = 2e19;         %edge density, if a parabolic profile is used
alpha_n         = 1;            %density profile parameter if a parabolic profile is used
lambda_n        = 0.05;         %tune exponential decay in edge (smaller = faster decay)
T0              = 5e3;          %core temp (eV)
T0_edge         = 1e2;          %edge temperature (eV)
alpha_T         = 1.5;          %profile shpae parameter for temperature if parabolic profiles are used
lambda_T        = 0.05;         %tune exponential decay in edge (smaller = faster decay)
iatanprofile    = false;        %use "advanced density profile"
fracap          = 0.01;         %tune "advanced density profile"
                   
%ions in plasma:   H       D          
n_species       = 2;              %number of ion species
A_ion           = [3       1      ];  %mass number of ions
Z_ion           = [2       1      ];  %charge of ions
conc_ion        = [0.12    0.88   ];  %concentration relative to ne
T_ion           = [5e3     5e3    ];   %ion temperature.
alpha_T_ion     = [alpha_T alpha_T];
alpha_n_ion     = [alpha_n alpha_n];


%ions in plasma:   H    D          He4
% n_species       = 3;              %number of ion species
% A_ion           = [1       2       4      ];  %mass number of ions
% Z_ion           = [1       1       2      ];  %charge of ions
% conc_ion        = [0.05    0.95    0.01  ];  %concentration relative to ne
% T_ion           = [5e3     5e3     1e6    ];   %ion temperature.
% alpha_T_ion     = [alpha_T alpha_T alpha_T];
% alpha_n_ion     = [alpha_n alpha_n alpha_n];

% n_species       = 7;            %number of ion species
% names           = [{'H'}   {'He3'} {'D'}   {'T'}   {'He4t'}  {'Be'}    {'He4'}];
% A_ion           = [1       3       2       3       4         9         4      ];  %mass number of ions
% Z_ion           = [1       2       1       1       2         4         2      ];  %charge of ions
% conc_ion        = [0.05    0.01    0.43    0.43    0.06      0.01      0.01   ];  %concentration relative to ne
% T_ion           = [T0      T0      T0      T0      T0        T0        400e3  ];  %ion temperature.
% alpha_T_ion     = [alpha_T alpha_T alpha_T alpha_T alpha_T   alpha_T   0.1    ];
% alpha_n_ion     = [alpha_n alpha_T alpha_n alpha_n alpha_T   alpha_n   1.0    ];





