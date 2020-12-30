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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  This is a small function to show how the product between the modified %
%  Bessel function of the first kind, and exp(-z) is approximated for    %
%  large z:  exp(-z)*besseli(a,z)                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
close all;
addpath('../simulator_functions'); % load, amongst others, bessel_times_exponent_approx.m

zr = 0:0.1:1000;    %real part of z(used in the plot)
zi = 0;           %imaginary part of z
z  = zr+1i*zi;      %total z
a = 5;

tic;
built_in_result = besseli(a,z).*exp(-z);
toc;
tic;
approximation   = bessel_times_exponent_approx(a,z);
toc;

%Plot the values obtained 
figure(1);
subplot(2,1,1)
semilogy(zr,abs(built_in_result));
xlim([min(zr) max(zr)]);
xlabel('Re(z)');
ylabel('I(a,z)*exp(-z)');
title('Result from MATLAB built-in functions');

subplot(2,1,2)
semilogy(zr,abs(approximation  ));
xlim([min(zr) max(zr)]);
xlabel('Re(z)');
ylabel('I(a,z)*exp(-z)');
title('Result from approximation');

figure(2);
semilogy(zr,abs(approximation-built_in_result));
xlim([min(zr) max(zr)]);
xlabel('Re(z)');
ylabel('error');
title('Absolute error');