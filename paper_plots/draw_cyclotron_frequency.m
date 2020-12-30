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
r0 = 2.97;       %major radius
a  = 0.95;       %minor radius
N  = 1001;       %number of gridpoints

mH = mp;
qH = e;

r = linspace(r0-a,r0+a,N).';
x = r-r0;
B = B0.*r0./(r0+x);
Omega_c = qH*B/mH;
f_c     = Omega_c/(2*pi);

%save background magnetic field
f1 = figure(1);
plot(x,B),hold on;
xlabel('$$x [\mathrm{m}]$$','interpreter','latex');
ylabel('$$B_z [\mathrm{T}]$$','interpreter','latex');
set(gca,'fontsize', 10)
savegoodplot('./plot_out/Bfield',[10 4],0)

%save hydrogen cyclotron frequency
f2 = figure(2);
plot(x,f_c/1e6),hold on;
plot(x,f*ones(N,1)/1e6,'-.');
legend('\Omega_{cH}','\omega');
xlabel('x [m]');
ylabel('f [MHz]');
title('Fundamental Hydrogen cyclotron harmonic and antenna frequency')
savegoodplot('./plot_out/cyclotronH',[10 4],0)


