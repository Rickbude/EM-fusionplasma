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

Nt = 20;
Np = 50;
Nr = 2;
R  = 2.97;
ap = 0.95;
a  = 1.15;

theta = linspace(0,2*pi,Nt);
phi   = linspace(0,pi,Np)';
r     = linspace(0,1,Nr)';

%wall torus
x1     = (R+a*cos(theta)).*cos(phi);
y1     = (R+a*cos(theta)).*sin(phi);
z1     = a*sin(theta).*ones(Np,1);

%plasma torus
x2     = (R+ap*cos(theta)).*cos(phi);
y2     = (R+ap*cos(theta)).*sin(phi);
z2     = ap*sin(theta).*ones(Np,1);

%left "lid"
x3     = -R+ap*cos(theta).*r;
y3     = zeros(Nr,Nt);
z3     = ap*sin(theta).*r;

%right "lid"
x4     =  R+ap*cos(theta).*r;
y4     = zeros(Nr,Nt);
z4     = ap*sin(theta).*r;

surf(x1,y1,z1,'facecolor',[166 166 166]/255), hold on
surf(x2,y2,z2,'facecolor',[247 163 195]/255), hold on
surf(x3,y3,z3,'linestyle','none','facecolor',[247 163 195]/255), hold on
surf(x4,y4,z4,'linestyle','none','facecolor',[247 163 195]/255), hold on
xlabel('r [m]');
zlabel('y [m]');

%draw 1D domain layout
rectangle('Position',[-4 -14  8 5],'FaceColor',[200 200 200]/255);    %rectanlge -a_s to a_s
rectangle('Position',[-3 -14  6 5],'linestyle','-.','linewidth',1.0,'FaceColor',[247 200 220]/255);    %rectanlge -a to a
rectangle('Position',[-2.5 -14  5 5],'FaceColor',[247 163 195]/255);    %rectanlge -a_p to a_p

%connecting lines between torus and 1D domain
plot3([R-a R+a R+a R-a],[0 0 0 0],[-0.2 -0.2 0.2 0.2],'linestyle','-.','linewidth',1,'color',[0 0 0]);    %rectanlge -a to a
plot3([R-a -3],[0 -14+5],[0 0],'-.','linewidth',1.0,'color',[0 0 0]);
plot3([R+a +3],[0 -14+5],[0 0],'-.','linewidth',1.0,'color',[0 0 0]);

%draw axis
quiver3(0,0,-a,0,0,3.5*a,'linewidth',1,'color',[0 0 0]);
quiver3(0,0,0,R+2*a,0,0 ,'linewidth',1,'color',[0 0 0]);
quiver3(-4.5,-14,0,10.5,0,0 ,'linewidth',1,'color',[0 0 0]);

%add text in 1D figure
text(-4.5,-15.2,'$-a_s$','interpreter','latex');
text( 3.9,-15.2, '$a_s$','interpreter','latex');
text(-3.5,-15.2,'$-a$','interpreter','latex');
text( 2.9,-15.2, '$a$','interpreter','latex');
text(-2.9,-15.2,'$-a_p$','interpreter','latex');
text( 2.4,-15.2, '$a_p$','interpreter','latex');
text(-0.1,-15.2,'$0$','interpreter','latex');

%add text to axis
text(4.8,-13,'$x$','interpreter','latex');
text(4.8,1 ,'$R$','interpreter','latex');

axis equal;
legend('metallic wall','plasma','Location','west','interpreter','latex');
set(gca, 'visible', 'off'); 
title('Simple Torus geometry');
grid on;
view([0 -1 0.2])
savegoodplot('./plot_out/1D_geom',[10 14],0)