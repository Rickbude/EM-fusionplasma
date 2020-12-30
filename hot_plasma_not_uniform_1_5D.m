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
addpath('hot_plasma_functions');         % load functions specifically for hot plasma (cp Stix tensor, etc)
addpath('cold_plasma_functions');        % load functions specifically for cold plasma (cp Stix tensor, etc)
addpath('EM_property_functions');        % load functions for EM properties (Poynting vector etc)
addpath('profile_functions');            % load functions for generating profiles
addpath('simulator_functions');          % load functions that support running simulations (general)
addpath('plasma_profiles');              % load parameter files
addpath('error_functions');              % load different error metrics

clear all;
close all;
profile on;
%load physical constants
constants;

tic;

%load "input_file.mat", which contains the name of the profile that should
%be loaded. If this file does not exist, or if you want to select a 
%different input file, run "load_input".
load('input_file');

%these should be moved to input files as well.
calc_power = 0;
compare_with_results_Dirk = 0;
compare_with_previous     = 0;
compare_with_all_orders   = 0;

%Bb: background magnetic field
[profiles,flags,settings] = load_profile(input_file);

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

run(input_file);

%% select below the solver you want to use
% hp_calc_E_truncated_fit: use truncated fit model
% hp_calc_E_truncated_taylor: use the truncated Taylor model
% hp_calc_E_all_orders: use the all-orders model
% uncomment the model of your choice.
[E,timings] = hp_calc_E_truncated_fit   (profiles,flags,settings);
%[E,timings] = hp_calc_E_truncated_taylor(profiles,flags,settings);
%[r_high_res,E,E_high_res,timings] = hp_calc_E_all_orders      (profiles,flags,settings);
%[E] = cp_calc_electric_field_1_5D(r,k,curve_k_parr,omega,Bb,m,n,q,Jant)
t_elapsed = toc;
fprintf('Done calculating fields (%.2f s), generating plots...\n',t_elapsed);

if(calc_power)
    P_abs  = zeros(N,n_species+1);
    w      = waitbar(0,sprintf('Calculating plasma current, species %i/%i',0,n_species+1));
    tprev  = toc;
    for species = 1:n_species+1
        J_p = hp_plasma_current_fast(r,E,k,curve_k_parr,omega,Bb,m(:,species),n(:,species),q(:,species),T(:,species),Jant);
        P_abs(:,species) = power_absorbed(J_p,E);
        waitbar(species/(n_species+1),w,sprintf('Calculating plasma current, species %i/%i',species,n_species+1));
    end
    close(w);

    figure(88)
    plot(r-r0,sum(P_abs,2)),hold on;
    plot(r-r0,P_abs),hold on;
    %set(gca,'YScale','log')
    xlabel('x [m]');
    ylabel('P_{abs} (a.u.)');
    legend('total','e','H','D');
    title('Power absorbtion based on the Poynting flux alone');

    figure(89)
    plot(r-r0,cumsum(sum(P_abs,2))),hold on;
    plot(r-r0,cumsum(P_abs));
    xlabel('r');
    ylabel('int P_{abs}');
    legend('total','e','H','D');
end

figure(20)


N_fft             = 2^(nextpow2(N)+3);
Fs                = N_fft/abs(max(r)-min(r))/(2*pi);    %samples per meter

spectrum        = abs(fftshift(fft(E,N_fft),1));
frequencies     = linspace(-Fs/2,Fs/2,N_fft)';

plot(frequencies,10*log10(spectrum./max(spectrum))),hold on;
plot([min(frequencies) max(frequencies)],[-30 -30],'--');
xlim([-2000 2000]);
ylim([-40 0]);
legend('E_x','E_y','E_z');
xlabel('k_x');
ylabel('log_{10}|E|');
title('spectrum of electric field components');

figure(21)

plot(frequencies,spectrum),hold on;
xlim([-2000 2000]);
legend('E_x','E_y','E_z');
xlabel('k_x');
ylabel('|E|');
title('spectrum of electric field components');

B = magnetic_field_1_5D(r,k,omega,E);

if(compare_with_results_Dirk)    
    load('results/results_dirk');
end

if(compare_with_previous && ~compare_with_all_orders)
    %load('results/previous_results');
    load('results/previous_results_dawson_exact');
end

if(compare_with_all_orders)
    %load('results/previous_results_all_orders_ITER_highres');
    %load('results/previous_results_truncated_taylor');
   % load('results/previous_results_ao_75H_5000');
   
   
   %%%%%   load('results/previous_results_ao_benchmark_5000');   %%%%%
     load('results/previous_results_ao_N9001_Nk9001_JET_75H');
    %load('results/He3_in_H_ao');
   
    %load('results/previous_results_all_orders_ITER');
    %load('results/previous_results_all_orders_no_Be_He4t');
    %interpolate previous results on current resolution
    E_interp = interp1(r_prev,E_prev,r,'spline');
    error    = E-E_interp;
end





%normalize Poyting flux, such that maximum = 1
if(iunityflux)
   [E,B] = normalize_flux(r,E,B);
end

Flux = Poynting_flux(r,E,B,curve_k_parr);



figure(1);

subplot(3,1,1);
plot(r,Bb);
xlabel('x [m]');
ylabel('B [T]');
title('Background magnetic field');
legend('B_x','B_y','B_z');


subplot(3,1,2);
plot(r,n),hold on;
xlabel('x [m]');
ylabel('n [m^{-3}]');
title('Density');
legend('n_e','n_H','n_D');

subplot(3,1,3);
plot(r,T*kB/e),hold on;
xlabel('x [m]');
ylabel('T [eV]');
title('Temperature');
legend('T_e','T_H','T_D');


figure(9)
mD          = m(:,3);
ne          = n(:,1);
omega_UH    = UH_frequency(Bb,ne);
omega_LH    = LH_frequency(Bb,ne,mD);
omega_Xmin  = Xmin_frequency (Bb,ne,mD);
omega_Xplus = Xplus_frequency(Bb,ne,mD);
omega_c     = cyclotron_frequency(q,Bb,m);

semilogy(r,abs(omega*ones(N,1))),hold on;
semilogy(r,omega_LH);
semilogy(r,omega_UH);
semilogy(r,omega_Xmin);
semilogy(r,omega_Xplus);
semilogy(r,real(omega_c));
title('Characteristic plasma frequencies');
xlabel('x [m]');
ylabel('\omega');


legend('\omega','\omega_{LH}','\omega_{UH}','\omega_{X_-}','\omega_{X_+}','\Omega_{ce}','\Omega_{cH}','\Omega_{cD}');


figure(10)

k_parr = k(:,3); %k_parr = kz
k_perp = cp_wavenumber_k_parr(omega,k_parr,Bb,m,n,q);

[~,ind]=min(abs(omega_c - real(omega)))
x = r - r0;


subplot(2,2,1);
root = real(k_perp(:,1).^2);
root_pos = zeros(N,1);
root_neg = zeros(N,1);
root_pos(root>0) =  root(root>0);
root_neg(root<0) = -root(root<0);
plot(x,root_pos,'b'),hold on;
plot(x,root_neg,'r'),hold on;
set(gca,'YScale','log')
title('Fast Wave');
xlabel('x [m]');
ylabel('Re(k_x^2)');
xline(x(ind(2)),'-.b');
legend('+','-','\Omega_{cH}');


subplot(2,2,2);
root = imag(k_perp(:,1).^2);
root_pos = zeros(N,1);
root_neg = zeros(N,1);
root_pos(root>0) =  root(root>0);
root_neg(root<0) = -root(root<0);
plot(r-r0,root_pos,'b'),hold on;
plot(r-r0,root_neg,'r'),hold on;
set(gca,'YScale','log')
title('Fast Wave');
xlabel('x [m]');
ylabel('Im(k_x^2)');
xline(x(ind(2)),'-.b');
legend('+','-','\Omega_{cH}');


subplot(2,2,3);
root = real(k_perp(:,2).^2);
root_pos = zeros(N,1);
root_neg = zeros(N,1);
root_pos(root>0) =  root(root>0);
root_neg(root<0) = -root(root<0);
plot(r-r0,root_pos,'b'),hold on;
plot(r-r0,root_neg,'r'),hold on;
set(gca,'YScale','log')
title('Slow Wave');
xlabel('x [m]');
ylabel('k');
ylabel('Re(k_x^2)');
xline(x(ind(2)),'-.b');
legend('+','-','\Omega_{cH}');



subplot(2,2,4);
root = imag(k_perp(:,2).^2);
root_pos = zeros(N,1);
root_neg = zeros(N,1);
root_pos(root>0) =  root(root>0);
root_neg(root<0) = -root(root<0);
plot(r-r0,root_pos,'b'),hold on;
plot(r-r0,root_neg,'r'),hold on;
set(gca,'YScale','log')
title('Slow Wave');
xlabel('x [m]');
ylabel('k');
ylabel('Im(k_x^2)');
xline(x(ind(2)),'-.b');
legend('+','-','\Omega_{cH}');

figure(55)

subplot(2,2,1);
root = real(k_perp(:,1).^2);
modified_logplot(x,root);
title('Fast Wave');
xlabel('x [m]');
ylabel('Re(k_x^2)');
xline(x(ind(2)),'-.b');
legend('k_x','\Omega_{cH}');


subplot(2,2,2);
root = imag(k_perp(:,1).^2);
modified_logplot(x,root);
title('Fast Wave');
xlabel('x [m]');
ylabel('Im(k_x^2)');
xline(x(ind(2)),'-.b');
legend('k_x','\Omega_{cH}');


subplot(2,2,3);
root = real(k_perp(:,2).^2);
modified_logplot(x,root);
title('Slow Wave');
xlabel('x [m]');
ylabel('k');
ylabel('Re(k_x^2)');
xline(x(ind(2)),'-.b');
legend('k_x','\Omega_{cH}');



subplot(2,2,4);
root = imag(k_perp(:,2).^2);
modified_logplot(x,root);
title('Slow Wave');
xlabel('x [m]');
ylabel('k');
ylabel('Im(k_x^2)');
xline(x(ind(2)),'-.b');
legend('k_x','\Omega_{cH}');



colors = get(gca,'colororder');

figure(2);
%subplot(2,1,1);

plot_components_and_abs(r-r0,real(E.*[1 1 100]),'-',colors);
if(compare_with_results_Dirk)
    plot_components_and_abs(r_Dirk-r0,real(E_Dirk).*[1 1 100],'--',colors);
end
if(compare_with_previous)
    plot_components_and_abs(r_prev-r0,real(E_prev).*[1 1 100],'--',colors);
end
%xlim([-a 0.8]);

plot_vertical_lines(ap,x_wall)    

h = flipud(findobj(gca,'Type','line'));

legend([h(1) h(2) h(3) h(4) h(6)],{'Re\{Ex\}','Re\{Ey\}','Re\{100*Ez\}','a_p','a'});
xlim([-x_wall x_wall])
xlabel('x [m]');
ylabel('E [V/m]')
%ylim([-20 20]);
%title(sprintf('D plasma, H minority, x_ant = %.2f, f = %.1fMHz, \\nu_\\omega =%.1g , n_D = %.0gm^{-3}, n_i = %.0gm^{-3}, B_0=%.2fT, k_z = %.1f',r0+x_a,f/1e6,nuom,n0,n0*conc_ion(1),B0,ntor/r0));
title('Real part of Electic field components for JET-parameters')

% 
% subplot(2,1,2);
% plot_components_and_abs(r-r0,imag(E.*[1 1 100]),'-',colors);
% if(compare_with_results_Dirk)
%     plot_components_and_abs(r_Dirk-r0,imag(E_Dirk).*[1 1 100],'--',colors);
% end
% 
% if(compare_with_previous)
%     plot_components_and_abs(r_prev-r0,imag(E_prev).*[1 1 100],'--',colors);
% end
% %xlim([-a 0.8]);
% legend('Im\{Ex\}','Im\{Ey\}','Im\{100*Ez\}');
% plot_vertical_lines(ap,x_wall)  
% xlabel('x [m]');
% ylabel('E [V/m]')
% %ylim([-20 20]);
% title(sprintf('D plasma, H minority, x_ant = %.2f, f = %.1fMHz, \\nu_\\omega =%.1g , n_D = %.0gm^{-3}, n_i = %.0gm^{-3}, B_0=%.2fT, k_z = %.1f',r0+x_a,f/1e6,nuom,n0,n0*conc_ion(1),B0,ntor/r0));


if(compare_with_previous)
    figure(5);
    subplot(1,3,1);
    plot(r,abs(E(:,1))),hold on;
    plot(r,abs(E_interp(:,1))),hold on;
    xlim(r0+[-x_wall x_wall]);
    xlim([-ap ap])
    xlabel('r[m]');
    ylabel('|E_x|');
    legend('modified linearized','all-orders');
    title('|E_x|, JET parameters');
    
    subplot(1,3,2);
    plot(r,abs(E(:,2))),hold on;
    plot(r,abs(E_interp(:,2))),hold on;
    xlim(r0+[-x_wall x_wall]);
    xlim([-ap ap])
    xlabel('r[m]');
    ylabel('|E_y|');
    legend('modified linearized','all-orders');
    title('|E_y|, JET parameters');
    
    subplot(1,3,3);
    plot(r,abs(E(:,3))),hold on;
    plot(r,abs(E_interp(:,3))),hold on;
    xlim(r0+[-x_wall x_wall]);
    xlim([-ap ap])
    xlabel('r[m]');
    ylabel('|E_z|');
    legend('modified linearized','all-orders');
    title('|E_z|, JET parameters');
   
end



figure(31)
subplot(3,1,1);
plot(r-r0,real(E(:,1))),hold on;
plot(r-r0,imag(E(:,1))),hold on;
xlabel('x[m]');
ylabel('E_x');
legend('real(E_x)','imag(E_x)');

subplot(3,1,2);
plot(r-r0,real(E(:,2))),hold on;
plot(r-r0,imag(E(:,2))),hold on;
xlabel('x[m]');
ylabel('E_y');
legend('real(E_y)','imag(E_y)');

subplot(3,1,3);
plot(r-r0,real(E(:,3))),hold on;
plot(r-r0,imag(E(:,3))),hold on;
xlabel('x[m]');
ylabel('E_z');
legend('real(E_z)','imag(E_z)');


% figure(21)
% subplot(2,1,1)
% mask = [1 1 100];
% for dim = 1:3
%    plot(r-r0,abs(E(:,dim))*mask(dim),'Color',colors(dim,:)),hold on;
% end
% title(sprintf('D plasma, H minority, x_ant = %.2f, f = %.1fMHz, \\nu_\\omega =%.1g , n_D = %.0gm^{-3}, n_i = %.0gm^{-3}, B_0=%.2fT, k_z = %.1f',r0+x_a,f/1e6,nuom,n0,n0*conc_ion(1),B0,ntor/r0));
% 
% 
% subplot(2,1,2)
% mask = [1 1 1];
% for dim = 1:3
%    plot(r-r0,angle(E(:,dim))*mask(dim),'Color',colors(dim,:)),hold on;
% end
% title(sprintf('D plasma, H minority, x_ant = %.2f, f = %.1fMHz, \\nu_\\omega =%.1g , n_D = %.0gm^{-3}, n_i = %.0gm^{-3}, B_0=%.2fT, k_z = %.1f',r0+x_a,f/1e6,nuom,n0,n0*conc_ion(1),B0,ntor/r0));




figure(3);
subplot(2,1,1);


plot_components_and_abs(r,real(omega)*real(B),'-',colors);
if(compare_with_results_Dirk)
    plot_components_and_abs(r_Dirk,real(omega)*real(B_Dirk),'--',colors);
end

if(compare_with_previous)
    plot_components_and_abs(r_prev,real(omega)*real(B_prev),'--',colors);
end

legend('Re\{\omega Bx\}','Re\{\omega By\}','Re\{\omega Bz\}');
xlabel('x [m]');
ylabel('B [A/m]')
%ylim([-20 20]);
title(sprintf('D plasma, H minority, x_ant = %.2f, f = %.1fMHz, \\nu_\\omega =%.1g , n_D = %.0gm^{-3}, n_i = %.0gm^{-3}, B_0=%.2fT, k_z = %.1f',r0+x_ant,f/1e6,nuom,n0,n0*conc_ion(1),B0,ntor/r0));

subplot(2,1,2);

plot_components_and_abs(r,real(omega)*imag(B),'-',colors);
if(compare_with_results_Dirk)
    plot_components_and_abs(r_Dirk,real(omega)*imag(B_Dirk),'--',colors);
end
if(compare_with_previous)
    plot_components_and_abs(r_prev,real(omega)*imag(B_prev),'--',colors);
end

legend('Im\{\omega Bx\}','Im\{\omega By\}','Im\{\omega Bz\}');
xlabel('x [m]');
ylabel('B [A/m]')
%ylim([-20 20]);
title(sprintf('D plasma, H minority, x_ant = %.2f, f = %.1fMHz, \\nu_\\omega =%.1g , n_D = %.0gm^{-3}, n_i = %.0gm^{-3}, B_0=%.2fT, k_z = %.1f',r0+x_ant,f/1e6,nuom,n0,n0*conc_ion(1),B0,ntor/r0));


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

%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(compare_with_all_orders)
    mask_zero   = abs((r-r0))>=ap;


    P_E_interp  = sum(abs(E_interp(~mask_zero,:)).^2);
    P_E         = sum(abs(E       (~mask_zero,:)).^2);
    P_error     = sum(abs(error   (~mask_zero,:)).^2);


    total_error = RRSE(E_interp(~mask_zero,:),E(~mask_zero,:));
    fprintf('total error:');
    disp(total_error)

    figure(140)
    semilogy(r-r0,abs(error)),hold on;
    %plot_vertical_lines(ap,x_wall);
    xline(-ap,'--')
    xline(ap,'--');
    legend('Error |Ex|','Error |Ey|','Error |Ez|');
    title('Absolute error');
    xlim([-x_wall x_wall]);
    xlabel('x[m]');
    ylabel('error');
    profile off;


end

figure(99);
x = r-r0;
plot(r,real(E(:,1)));
xlim(r0+[-x_wall x_wall]);
title('Electric field for ICRH in JET')
xlabel('Distance to device center [m]')
ylabel('Real part of Ex [V/m]');
[~,ind]=min(abs(omega_c - omega))
xline(r(ind(2)),'-.b');
legend('Electric field','Absorbtion layer');

figure(100);
x = r-r0;
plot(r,real(E(:,1)));
xlim([r0-x_wall r0]);
title('Electric field with the slow model')
xlabel('Distance to device center [m]')
ylabel('Real part of Ex [V/m]');
[~,ind]=min(abs(omega_c - omega))
xline(x(ind(2)),'-.b');

figure(101);
x = r-r0;
plot(r,real(E(:,1)));
xlim([r0-x_wall r0+x_wall]);
title('Simulated Electric Field')
xlabel('Distance to device center [m]')
ylabel('Real part of Ex [V/m]');



E_prev = E;
B_prev = B;
r_prev = r;
save('results/previous_results','r_prev','E_prev','B_prev');
