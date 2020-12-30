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
addpath('../plasma_profiles');              % load parameter files

clearvars;
close all;
profile on;
%load physical constants
constants;
profile_name = 'JET_test_Dirk_75H';

N = 2001;
N_orders = 12;
x_test = -0.10;
r_test   = 2.97+x_test;

%Bb: background magnetic field
[profiles,flags,settings] = load_profile(profile_name);

Bb = profiles.Bb;
n  = profiles.n;
m  = profiles.m;
q  = profiles.q;
T  = profiles.T;
r  = profiles.r;
k  = profiles.k;
omega = settings.omega;
r0 = mean(r);


[~,r_test_ind] = min(abs(r-r_test));

Bb = Bb(r_test_ind,:);
n  = n (r_test_ind,:);
m  = m (r_test_ind,:);
q  = q (r_test_ind,:);
T  = T (r_test_ind,:);

kx    = linspace(-1000,1000,N).';
ky    = ones(N,1)*mean(k(:,2));
kz    = ones(N,1)*mean(k(:,3));
kz0   = mean(kz);

rho_L = larmour_radius(T,m,q,Bb);
rho_L = rho_L(2);

k     = [kx ky kz];

dk    = 1;
h     = (max(kx)-min(kx))/(N-1);
k_test= [0 0 kz0];
kxt   = k_test(1);
dkx   = dk*[1 0 0];
dkz   = dk*[0 0 1];

sigmas   = hp_dielectric_tensor_swanson(N_orders,omega,k,Bb,T,q,n,m);

[epsilon, d_epsilon, H_epsilon] = hp_dielectric_tensor_and_derivatives(N_orders,omega,k_test.*ones(N,1),dk,Bb,T,q,n,m);


k_shift  = k-k_test;

scriptE0 = epsilon - kxt.*squeeze(d_epsilon(:,:,:,1)) + 1/2.*kxt.^2.*squeeze(H_epsilon(:,:,:,1,1)); 
scriptE1 = squeeze(d_epsilon(:,:,:,1)) - kxt.*squeeze(H_epsilon(:,:,:,1,1));
scriptE2 = 0.5*squeeze(H_epsilon(:,:,:,1,1));

epsilon_approx = scriptE0 + scriptE1.*kx + scriptE2.*kx.^2;

figure(2)
%title(sprintf('dielectric tensor, linearized around k = %i, N',kz0));
colors = get(gca,'colororder');

plot(kx*rho_L,real(sigmas(:,1,1)),'Color',colors(1,:)), hold on;
plot(kx*rho_L,real(epsilon_approx(:,1,1)),'--','Color',colors(1,:)), hold on;

xlabel('$k_x \rho_{\mathrm{LH}} [-]$','interpreter','latex')
ylabel('$\Re(\epsilon_{11})$','interpreter','latex');
xlim([ -3 3]);
ylim([ -300 2000]);
legend('calculated','2nd order Taylor','location','southeast');

%title(sprintf('Approximation of \\epsilon_{11} at x = %.2f',r_test-r0));
 
savegoodplot('./plot_out/eps11_taylor',[10 6], 0)

profile off;
