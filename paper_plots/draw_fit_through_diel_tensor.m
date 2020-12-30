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

%load libraries
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

Nk       = 100;     %number of samples in 'k-space' for the fit
N_orders = 12;      %maximum harmonic to include in hot plasma tensor 
x_test   = -0.1;    %sample location of the dielectric tensor
kx_limit = 1000;    %width of the spectrum to conside
Np       = 6;       %polynomial order

%Load profile properties
%Bb: background magnetic field.
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


%find closest matching gridpoint
r0              = (max(r)+min(r))/2;
r_test          = x_test+r0;
[~,r_test_ind]  = min(abs(r-r_test));

%take background properties on closest matching gridpoint
Bb = Bb(r_test_ind,:);
n  = n (r_test_ind,:);
m  = m (r_test_ind,:);
q  = q (r_test_ind,:);
T  = T (r_test_ind,:);

kx_compressed = linspace(-1,1,Nk)';
kx    = kx_compressed*kx_limit;
ky    = ones(Nk,1)*k(r_test_ind,2);
kz    = ones(Nk,1)*k(r_test_ind,3);
kz0   = k(r_test_ind,3);
k     = [kx ky kz];




epsilon= zeros(Nk,3,3);
for ind = 1:Nk
    k_t = k(ind,:);
    epsilon(ind,:,:) = squeeze(hp_dielectric_tensor_swanson(N_orders,omega,k_t,Bb,T,q,n,m));
end

rho_L = larmour_radius(T,m,q,Bb);
rho_L = rho_L(2);
ak = zeros(3,3,Np+1);

for dim1 = 1:3
    for dim2 = 1:3
        p = polyfit(kx_compressed(abs(kx)<1/rho_L),epsilon(abs(kx)<1/rho_L,dim1,dim2),Np);     
        ak(dim1,dim2,:) = p./(kx_limit.^fliplr(0:Np));
    end
end




figure(88)
colors = get(gca,'colororder');
plot(kx*rho_L,squeeze(real(epsilon(:,1,1))),'Color',colors(1,:)),hold on;


p = polyfit(rho_L*kx(kx>0),epsilon(kx>0,1,1),6);
actual = epsilon(kx>0,1,1);
fitted = polyval(squeeze(ak(1,1,:)),kx);
plot(kx*rho_L,real(fitted),'--','Color',colors(1,:));
plot([1 1],[min(real(actual)) max(real(actual))*2],'-.','Color',[0 0 0])
plot(-[1 1],[min(real(actual)) max(real(actual))*2],'-.','Color',[0 0 0])
ylim([min(real(actual)) max(real(actual))]);

% xlabel('k_x \rho_{LH} [-]');
% ylabel('Re(\epsilon_{11})');
% title(sprintf('\\epsilon_{11} at x = %.2f m',x_test));
% legend('actual','6th order fit','fit bounds');

xlabel('$k_x \rho_{\mathrm{LH}} [-]$','interpreter','latex')
ylabel('$\Re(\epsilon_{11})$','interpreter','latex');
xlim([ -3 3]);
ylim([ -300 2000]);
legend('actual','6th order fit','fit bounds','location','southeast');

savegoodplot('./plot_out/eps11_fit',[10 6],0)


figure(1)

colors = get(gca,'colororder');
plotnum = 1;
counter = 1;
for dim1 = 1:3
    for dim2 =1:3
        subplot(3,3,plotnum);
        
        plot(rho_L*kx,real(epsilon(:,dim1,dim2)),'Color',colors(counter,:)), hold on;
        
        
        p = polyfit(rho_L*kx(kx>0),epsilon(kx>0,dim1,dim2),6);
        plot(rho_L*kx,real(polyval(squeeze(ak(dim1,dim2,:)),kx)),'-.r');
        
       
        counter = counter+1;
        if(counter>length(colors))
            counter = 1;
        end
        plotnum=plotnum+1;
        xlabel('\rho_L k_x [-]')
        ylabel(sprintf('Re(\\epsilon_{%i%i})',dim1,dim2));
        ylim([min(real(epsilon(:,dim1,dim2)))-1 max(real(epsilon(:,dim1,dim2)))+1]);
        %legend('original','fit');
    end
end
subplot(3,3,1)
profile off;

