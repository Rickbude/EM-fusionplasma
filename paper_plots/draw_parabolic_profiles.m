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

clearvars;
close all;

f  = 51e6;       %antenna frequency
e  = 1.602e-19;  %elementary charge
mp = 1.673e-27;  %proton mass
B0 = 3.45;       %magnetic field on axis
r0 = 2.97;
a  = 1.05;
ap = 0.95;
N  = 1001;

mH = mp;
qH = e;

r = linspace(r0-a,r0+a,N).';
x = r-r0;
B = B0.*r0./(r0+x);
Omega_c = qH*B/mH;
f_c     = Omega_c/(2*pi);

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

n = parabolic_profile(x,ap,n0,n0_edge,alpha_n,lambda_n);
T = parabolic_profile(x,ap,T0,T0_edge,alpha_T,lambda_T);

gcf = figure(1);
set(gcf,'position',[200 200 400 300])
subplot(2,1,1)
plot(x,n/1e19)
xlabel('$$x[\mathrm{m}]$$','interpreter','latex')
ylabel('$$n_e[10^{19} \mathrm{m^{-3}}]$$','interpreter','latex');
%title('Electron density');
xline(-ap,'-.');
xline(ap,'-.');
legend('$n_e$','$a_p$','Interpreter','latex')

set(gca,'FontSize', 10);
set(gca,'fontname','times');

grid on;

subplot(2,1,2);
plot(x,T/1e3);
xlabel('$$x[\mathrm{m}]$$','interpreter','latex')
ylabel('$$T_e[\mathrm{keV}]$$','interpreter','latex');
%title('Electron temperature');
xline(-ap,'-.');
xline(ap,'-.');
legend('$T_e$','$a_p$','Interpreter','latex');
set(gca,'FontSize', 10);
set(gca,'fontname','times');

% set(gca,'fontname', 'times')

savegoodplot('./plot_out/T_and_n',[10 8],0)


saveas(gcf,'test.png');