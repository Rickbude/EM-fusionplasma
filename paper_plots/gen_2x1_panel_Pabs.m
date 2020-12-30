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

addpath('../plasma_property_functions');    % load functions that contain  plasma properties (cyclotron frequency etc)
addpath('../hot_plasma_functions_clean');         % load functions specifically for hot plasma (cp Stix tensor, etc)
addpath('../cold_plasma_functions');        % load functions specifically for cold plasma (cp Stix tensor, etc)
addpath('../EM_property_functions');        % load functions for EM properties (Poynting vector etc)
addpath('../profile_functions');            % load functions for generating profiles
addpath('../simulator_functions');          % load functions that support running simulations (general)
addpath('../error_functions');              % load different error metrics

clear all;
close all;
profile on;
%load physical constants
constants;

tic;
x_min = -1.05;
x_max = 0;
dim   = 1;

load('../input_file');
%profile_name = 'JET_test_Dirk_75H';
input_file = strcat('../',input_file);
[profiles,flags,settings] = load_profile(input_file);
%Bb: background magnetic field

%calculate the electric field by using the truncated Taylor model
E_taylor    = hp_calc_E_truncated_taylor   (profiles,flags,settings);
%calculate the electric field by using the all-orders model
[~,E_ao,~]   = hp_calc_E_all_orders      (profiles,flags,settings);

Bb           = profiles.Bb;
k            = profiles.k;
m            = profiles.m;
n            = profiles.n;
q            = profiles.q;
T            = profiles.T;
Jant         = profiles.Jant;
r            = profiles.r;
omega        = settings.omega;
curve_k_parr = flags.curve_k_parr;
n_species    = size(n,2);


disp('Calculating power');
Pabs_ao      = zeros(settings.N,n_species);
Pabs_taylor  = zeros(settings.N,n_species);
tprev  = toc;
for species = 1:n_species
    fprintf('Calculating species %i out of %i\n',species,n_species);
    
    m_j = m(:,species);
    n_j = n(:,species);
    q_j = q(:,species);
    T_j = T(:,species);
    
    Jp_taylor = hp_plasma_current_fast(r,E_taylor,k,curve_k_parr,omega,Bb,m_j,n_j,q_j,T_j,Jant);
    Jp_ao     = hp_plasma_current_fast(r,E_ao    ,k,curve_k_parr,omega,Bb,m_j,n_j,q_j,T_j,Jant);
    Pabs_taylor(:,species) = power_absorbed(Jp_taylor,E_taylor);
    Pabs_ao (:,species)    = power_absorbed(Jp_ao    ,E_ao    );
end



r0 = mean(profiles.r);
x  = profiles.r - r0;

ap = 0.95;
x_wall = settings.x_wall;
x_bnd  = bitand(x>x_min,x<x_max);

figure(88)

subplot(2,1,1)
plot(r-r0,Pabs_ao),hold on;
plot(r-r0,sum(Pabs_ao,2),'-.'),hold on;
set(gca,'FontSize', 10);
set(gca,'fontname','times');
title('Absorbed power, all-orders model','interpreter','latex');
xlabel('$x [\mathrm{m}]$','interpreter','latex');
ylabel('$P_{abs} [\mathrm{a.u.}]$','interpreter','latex');
legend('e','H','D','total');
grid on;

subplot(2,1,2)
plot(r-r0,Pabs_taylor),hold on;
plot(r-r0,sum(Pabs_taylor,2),'-.'),hold on;
set(gca,'FontSize', 10);
set(gca,'fontname','times');
title('Absorbed power, truncated taylor model','interpreter','latex');
xlabel('$x [\mathrm{m}]$','interpreter','latex');
ylabel('$P_{abs} [\mathrm{a.u.}]$','interpreter','latex');
legend('e','H','D','total');
grid on;

savegoodplot('./plot_out/compare_Pabs',[10 12],0)