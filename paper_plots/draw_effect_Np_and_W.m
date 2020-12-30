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

load('../input_file');
%profile_name = 'JET_test_Dirk_75H';
input_file = strcat('../',input_file);
[profiles,flags,settings] = load_profile(input_file);
%Bb: background magnetic field

all_W  = linspace(1e-3,4,200);
all_Np = 2:2:12;

%load high-resolution results, and interpolate to the lower resolution
load('results/previous_results_ao_N12601_Nk12601_JET_75H');
E_interp  = interp1(r_prev,E_prev,profiles.r,'spline');

%Only calculate the error for the entries within the plasma
inside_plasma = abs((profiles.r-settings.r0))<=settings.ap;

total_errors = zeros(length(all_Np),length(all_W));
total_times  = zeros(length(all_Np),length(all_W));

indNp = 1;
for Np = all_Np 
    indW  = 1;
    for W = all_W
        %reload the profiles with overloaded settings
        override.zeta = W;
        override.Np   = Np;
        override.Nk   = 50;
        [profiles,flags,settings] = load_profile_override(input_file,override);
        
        %reset the timer to 0, calculate E-field, store time
        tic;
        [E,timings] = hp_calc_E_truncated_fit   (profiles,flags,settings);
        total_times( indNp,indW) = timings.end;
        
        %calculate and store RRSE        
        total_error = RRSE(E_interp(inside_plasma,:),E(inside_plasma,:));
        total_errors(indNp,indW) = mean(total_error);
        
        indW = indW+1;        
    end
    indNp = indNp+1;
end

E_taylor     = hp_calc_E_truncated_taylor(profiles,flags,settings);
error_taylor = mean(RRSE(E_interp(inside_plasma,:),E_taylor(inside_plasma,:)));

override.Nk   = 200;
[profiles,flags,settings] = load_profile_override(input_file,override);
[~,E_ao,~]   = hp_calc_E_all_orders      (profiles,flags,settings);
error_ao     = mean(RRSE(E_interp(inside_plasma,:),E_ao(inside_plasma,:)));

%Make a plot of error as function of W, for different polynomial orders
figure(1);
semilogy(all_W,total_errors),hold on;
yline(error_taylor,'-.');
yline(error_ao,'-.');
hold off;
ylim([1e-3 1]);
legend('$N_p=2$','$N_p=4$','$N_p=6$','$N_p=8$','$N_p=10$','$N_p=12$','Location','southeast','interpreter','latex');
xlabel('Fit window $W$','interpreter','latex');
ylabel('RRSE','interpreter','latex');
savegoodplot('./plot_out/effect_Np_and_W',[10 6],0);

%Make a plot of time taken as function of W, for different polynomial
%orders
figure(2);
plot(all_W,total_times);
legend('N_p=2','N_p=4','N_p=6','N_p=8','N_p=10','N_p=12','Location','southeast');
xlabel('Fit window W');
ylabel('Time taken [s]');
savegoodplot('./plot_out/time_Np_and_W',[10 6],0);