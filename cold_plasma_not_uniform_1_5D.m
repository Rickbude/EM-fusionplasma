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

addpath('plasma_property_functions');    % load functions that contain  plasma properties (cyclotron frequency etc)
addpath('cold_plasma_functions');        % load functions for cold plasma (cp Stix tensor, etc)
addpath('EM_property_functions');        % load functions for EM properties (Poynting vector etc)
addpath('profile_functions');
addpath('plasma_profiles');
clearvars;
close all;
profile on;
%load "constants of nature"
constants;

tic;
theta   = pi/2;   

%load('results_dirk');

%profile_name = 'JET';

%Bb: background magnetic field
[omega,r,Bb,n,m,q,Jant,k] = load_profile(profile_name);

run(strcat('plasma_profiles/',profile_name,'.m'));

E = cp_calc_electric_field_1_5D(r,k,curve_k_parr,omega,Bb,m,n,q,Jant);


B = magnetic_field_1_5D(r,k,omega,E);

Flux = Poynting_flux(r,E,B);

t_elapsed = toc;
fprintf('Done calculating (%.0f ms), generating plots...\n',t_elapsed*1000);

%normalize Poyting flux, such that maximum = 1
if(false)
    dFluxx = max(real(Flux(:,1)))-min(real(Flux(:,1)));
    E = E./sqrt(dFluxx);
    B = B./sqrt(dFluxx);
    Flux = Poynting_flux(r,E,B);
end

figure(1);

subplot(2,1,1);
plot(r,Bb);
xlabel('x [m]');
ylabel('B [T]');
title('Background magnetic field');
legend('B_x','B_y','B_z');


subplot(2,1,2);
plot(r,n),hold on;
%plot(r_Dirk,n_Dirk);
xlabel('x [m]');
ylabel('n [m^{-3}]');
title('density');
legend('ne','n_H','n_D');


figure(9)
subplot(2,1,1);
mi          = m(:,3);
ne          = n(:,1);
ni          = n(:,3);
omega_UH    = UH_frequency(Bb,ne);
omega_LH    = LH_frequency(Bb,ne,mi);
omega_Xmin  = Xmin_frequency (Bb,ne,mi);
omega_Xplus = Xplus_frequency(Bb,ne,mi);
omega_c     = cyclotron_frequency(q,Bb,m);
plot(r,log10(abs(omega*ones(N,1)))),hold on;
plot(r,log10(omega_LH));
plot(r,log10(omega_UH));
plot(r,log10(omega_Xmin));
plot(r,log10(omega_Xplus));
plot(r,log10(real(omega_c)));
title('Characteristic plasma frequencies');
xlabel('x [m]');
ylabel('\omega');


legend('\omega','\omega_{LH}','\omega_{UH}','\omega_{X_-}','\omega_{X_+}','\Omega_{ce}','\Omega_{cH}','\Omega_{cD}');

subplot(2,1,2);
k      = cp_wavenumber(omega,theta,Bb,m,n,q);
plot(r,real(k.^2));
title('wavenumber');
xlabel('x [m]');
ylabel('k');

load('results_dirk');
colors = get(gca,'colororder');

figure(2);
subplot(2,1,1);

plot_components_and_abs(r,real(E),'-',colors);
plot_components_and_abs(r_Dirk,real(E_Dirk),'--',colors);

legend('Re\{Ex\}','Re\{Ey\}','Re\{Ez\}','|E|');
xlabel('x [m]');
ylabel('E [V/m]')
%ylim([-20 20]);
title(sprintf('D plasma, H minority, x_a = %.2f, f = %.1fMHz, \\nu_\\omega =%.1g , n_D = %.0gm^{-3}, n_i = %.0gm^{-3}, B_0=%.2fT, k_z = %.1f',r0+x_a,f/1e6,nuom,n0,n0*conc_ion(1),B0,ntor/r0));


subplot(2,1,2);
plot_components_and_abs(r,imag(E),'-',colors);
plot_components_and_abs(r_Dirk,imag(E_Dirk),'--',colors);


legend('Im\{Ex\}','Im\{Ey\}','Im\{Ez\}','|E|');
xlabel('x [m]');
ylabel('E [V/m]')
%ylim([-20 20]);
title(sprintf('D plasma, H minority, x_a = %.2f, f = %.1fMHz, \\nu_\\omega =%.1g , n_D = %.0gm^{-3}, n_i = %.0gm^{-3}, B_0=%.2fT, k_z = %.1f',r0+x_a,f/1e6,nuom,n0,n0*conc_ion(1),B0,ntor/r0));




figure(3);
subplot(2,1,1);


plot_components_and_abs(r,real(omega)*real(B),'-',colors);
plot_components_and_abs(r_Dirk,real(omega)*real(B_Dirk),'--',colors);


legend('Re\{\omega Bx\}','Re\{\omega By\}','Re\{\omega Bz\}','|\omega B|');
xlabel('x [m]');
ylabel('B [A/m]')
%ylim([-20 20]);
title(sprintf('D plasma, H minority, x_a = %.2f, f = %.1fMHz, \\nu_\\omega =%.1g , n_D = %.0gm^{-3}, n_i = %.0gm^{-3}, B_0=%.2fT, k_z = %.1f',r0+x_a,f/1e6,nuom,n0,n0*conc_ion(1),B0,ntor/r0));

subplot(2,1,2);

plot_components_and_abs(r,real(omega)*imag(B),'-',colors);
plot_components_and_abs(r_Dirk,real(omega)*imag(B_Dirk),'--',colors);

legend('Im\{\omega Bx\}','Im\{\omega By\}','Im\{\omega Bz\}','\omega|B|');
xlabel('x [m]');
ylabel('B [A/m]')
%ylim([-20 20]);
title(sprintf('D plasma, H minority, x_a = %.2f, f = %.1fMHz, \\nu_\\omega =%.1g , n_D = %.0gm^{-3}, n_i = %.0gm^{-3}, B_0=%.2fT, k_z = %.1f',r0+x_a,f/1e6,nuom,n0,n0*conc_ion(1),B0,ntor/r0));


figure(4)
subplot(2,1,1);
plot(r,real(Flux(:,1)));
xlabel('r [m]');
ylabel('S_p [W/m]');
title('Real Poyinting flux in -x direction');

subplot(2,1,2);
plot(r,imag(Flux(:,1)));
xlabel('r [m]');
ylabel('S_p [W/m]');
title('Imagionary Poyinting flux in -x direction');

profile off;
