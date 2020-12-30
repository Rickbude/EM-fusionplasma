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
addpath('../hot_plasma_functions');         % load functions specifically for hot plasma (cp Stix tensor, etc)
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

dim   = 1;

load('../input_file');
%profile_name = 'JET_test_Dirk_75H';
input_file = strcat('../',input_file);




%Bb: background magnetic field

override.zeta = 0.01;
override.Np   = 2;
[profiles,flags,settings] = load_profile_override(input_file,override);
E_fit_low = hp_calc_E_truncated_fit   (profiles,flags,settings);

override.zeta = 1.2;
override.Np   = 8;
[profiles,flags,settings] = load_profile_override(input_file,override);
E_fit_high = hp_calc_E_truncated_fit   (profiles,flags,settings);


E_taylor   = hp_calc_E_truncated_taylor(profiles,flags,settings);
[~,E_ao,~] = hp_calc_E_all_orders      (profiles,flags,settings);

x  = profiles.r - mean(profiles.r);
ap = 0.95;
x_wall = settings.x_wall;
x_min = -settings.ap;
x_max =  0;
x_bnd  = bitand(x>x_min,x<x_max);

load('results/previous_results_ao_N12601_Nk12601_JET_75H');
E_interp  = interp1(r_prev,E_prev,profiles.r,'spline');

figure(1);
colors = get(gca,'colororder');
subplot(3,1,1);

plot(x(x_bnd),abs(E_ao(x_bnd,dim)));
xline(-ap,'--');   
xlim([x_min x_max])
set(gca,'FontSize', 10);
set(gca,'fontname','times');
xlabel('$x [\mathrm{m}]$','interpreter','latex');
ylabel('$|E_x| [\mathrm{V/m}]$','interpreter','latex')
error = mean(RRSE(E_interp(x_bnd,dim),E_ao(x_bnd,dim)));
title(sprintf('All-orders, RRSE = %.2e',error),'interpreter','latex');
grid on

subplot(3,1,2);
plot(x(x_bnd),abs(E_taylor(x_bnd,dim)));
xline(-ap,'--');   
xlim([x_min x_max])
set(gca,'FontSize', 10);
set(gca,'fontname','times');
xlabel('$x [\mathrm{m}]$','interpreter','latex');
ylabel('$|E_x| [\mathrm{V/m}]$','interpreter','latex')
error = mean(RRSE(E_interp(x_bnd,dim),E_taylor(x_bnd,dim)));
title(sprintf('Truncated Taylor: RRSE = %.2e',error),'interpreter','latex');
grid on


subplot(3,1,3);
plot(x(x_bnd),abs(E_fit_high(x_bnd,dim)));
xline(-ap,'--');   
xlim([x_min x_max])
set(gca,'FontSize', 10);
set(gca,'fontname','times');
xlabel('$x [\mathrm{m}]$','interpreter','latex');
ylabel('$|E_x| [\mathrm{V/m}]$','interpreter','latex')
error = mean(RRSE(E_interp(x_bnd,dim),E_fit_high(x_bnd,dim)));
title(sprintf('Truncated polynomial: RRSE = %.2e',error),'interpreter','latex');
grid on

savegoodplot('./plot_out/Ex_compare_RRSE',[10 14],0)



