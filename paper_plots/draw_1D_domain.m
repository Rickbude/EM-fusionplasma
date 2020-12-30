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

ap = 0.95;
a  = 1.05;
as = 1.35;
d  = 0.1;

figure(1);
hold on;

%draw 1D domain layout
rectangle('Position',[-as 0  2*as 0.2],'FaceColor',[200 200 200]/255);    %rectanlge -a_s to a_s
rectangle('Position',[-a  0  2*a  0.2],'FaceColor',[247 200 220]/255);    %rectanlge -a to a
rectangle('Position',[-ap 0  2*ap 0.2],'FaceColor',[247 163 195]/255);    %rectanlge -a_p to a_p

%add markers
plot([-as -as],[-0.01 0.01],'linewidth',1,'color',[0 0 0]);
plot([-a  -a ],[-0.01 0.01],'linewidth',1,'color',[0 0 0]);
plot([-ap -ap],[-0.01 0.01],'linewidth',1,'color',[0 0 0]);
plot([ 0   0 ],[-0.01 0.01],'linewidth',1,'color',[0 0 0]);
plot([ ap  ap],[-0.01 0.01],'linewidth',1,'color',[0 0 0]);
plot([ a   a ],[-0.01 0.01],'linewidth',1,'color',[0 0 0]);
plot([ as  as],[-0.01 0.01],'linewidth',1,'color',[0 0 0]);

%add marker text
text(-as-0.20 ,-0.05,'$-a_s$','interpreter','latex');
text( as-0.05 ,-0.05, '$a_s$','interpreter','latex');
text(-a -0.15 ,-0.05,'$-a$'  ,'interpreter','latex');
text( a       ,-0.05, '$a$'  ,'interpreter','latex');
text(-ap-0.05 ,-0.05,'$-a_p$','interpreter','latex');
text( ap-0.05 ,-0.05, '$a_p$','interpreter','latex');
text(-0.03     ,-0.05,'$0$'   ,'interpreter','latex');

%add "plasma text"
text(-0.3,0.1,'Core Plasma','interpreter','latex');
text(-as+0.01,0.13,'Spill-','interpreter','latex');
text(-as+0.01,0.07,'over'  ,'interpreter','latex');
text(a+0.01,0.13,'Spill-','interpreter','latex');
text(a+0.01,0.07,'over'  ,'interpreter','latex');

%add axis
annotation('arrow',[0.12 0.95],[0.32 0.32]);
text(as+0.1,0.05,'$x$','interpreter','latex');

xlim([-as-0.05 as+0.05]);
ylim([-0.1 0.35]);

set(gca,'FontSize', 10);
set(gca,'fontname','times');

axis off;

savegoodplot('./plot_out/1D_geom',[10 4],0)

