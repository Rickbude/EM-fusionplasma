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
%  This is a small function to show how well the plasma dispersion       %
%  function is approximated. You will see that the approximations are    %
%  best for real values of zeta; and the results can diverge quickly     %
%  when the imaginary part of zeta becomes larger. Unfortunately, this   %
%  is also the region where the built-in dawson function is slowest.     %
%  For these regions, an additional expansion should be used.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
close all;
addpath('../simulator_functions'); % load, amongst others, bessel_times_exponent_approx.m

N  = 1000;              %samples
t  = linspace(0,1,N)';  %a "time"
zr = t*20;              %real part of z(used in the plot)
zi = 0.1;               %imaginary part of z
Nf = 5;                 %Terms in the continued fraction
z  = zr+1i*zi;          %total z

% The plasma dispersion function can be related to the error function, see 
% Fried & Conte 1961: the plasma dispersion function:
% Z = i*sqrt(pi) * exp(-zeta^2) * (1 + erf(i*z))
% Calculate this using built-in functions
tic;
built_in_result_1 = 1i*sqrt(pi).*exp(-z.^2).*(1+1i*erfi(z));
toc;
%It may also be written in terms of the Dawson function, which is slightly
%better behaving numerically than the erfi function
tic;
built_in_result_2 = 1i.*sqrt(pi).*exp(-z.^2) - 2*dawson(z);
toc;

tic;
%Now calculate it using the set of approximations used in the paper 
approximation_1   = plasma_dispersion_function(z);
toc;
tic;
%Calculate the dispersion function using the continued fractions from Fried
%and Conte. This works especially well for Im(zeta)>1, for somewhat larger
%Re(zeta)
approximation_2   = plasma_dispersion_function_continued_fraction(z,Nf);
toc;

%Plot the values obtained 
figure(1);
subplot(2,2,1)
plot(t,(real(built_in_result_2))),hold on;
plot(t,(real(built_in_result_1))),hold off;
xlim([min(t) max(t)]);
xlabel('t');
ylabel('|real(Z(zeta(t)))|');
legend('dawson','erfi');
title('Result from MATLAB built-in functions');

subplot(2,2,2)
plot(t,(imag(built_in_result_2))),hold on;
plot(t,(imag(built_in_result_1))),hold off;
xlim([min(t) max(t)]);
xlabel('t');
ylabel('|imag(Z(zeta(t)))|');
legend('dawson','erfi');
title('Result from MATLAB built-in functions');

subplot(2,2,3)
plot(t,(real(approximation_1))),hold on;
plot(t,(real(approximation_2))),hold off;
xlim([min(t) max(t)]);
xlabel('t');
ylabel('|real(Z(zeta(t)))|');
title('Result from approximation');
legend('composite','continued fraction');

subplot(2,2,4)
plot(t,(imag(approximation_1))),hold on;
plot(t,(imag(approximation_2))),hold on;
xlim([min(t) max(t)]);
xlabel('t');
ylabel('|imag(Z(zeta(t)))|');
title('Result from approximation');
legend('composite','continued fraction');

figure(2);
semilogy(zr,abs(built_in_result_1-built_in_result_2));
xlim([min(t) max(t)]);
xlabel('t');
ylabel('error');
title('Difference between built-in functions');

figure(3);
semilogy(t,abs(approximation_1-built_in_result_2));
xlim([min(t) max(t)]);
xlabel('Re(z)');
ylabel('error');
title('Absolute error');

figure(4);
semilogy(zr,abs(real(z))),hold on;
semilogy(zr,abs(imag(z))),hold on;
legend('abs(real(zeta(t)))','abs(imag(zeta(t)))');
title('zeta');
xlabel('t');
ylabel('zeta');
