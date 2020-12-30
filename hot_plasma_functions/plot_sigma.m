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

function [] = plot_sigma(r,sigmas,sigmas_dx,sigmas_dxx)
%PLOT_SIGMA Summary of this function goes here
%   Detailed explanation goes here
    
    hold off;
    figure();
    plot(r,sigmas(:,1,1),'Color',[1,0  ,0]),hold on;
    plot(r,sigmas(:,1,2),'Color',[1,0.3,0]),hold on;
    plot(r,sigmas(:,1,3),'Color',[1,0.6,0]),hold on;

    plot(r,sigmas(:,2,1),'Color',[0.6,0.8,0  ]),hold on;
    plot(r,sigmas(:,2,2),'Color',[0  ,0.8,0  ]),hold on;
    plot(r,sigmas(:,2,3),'Color',[0  ,0.8,0.6]),hold on;

    plot(r,sigmas(:,3,1),'Color',[0  ,0,1]),hold on;
    plot(r,sigmas(:,3,2),'Color',[0.3,0,1]),hold on;
    plot(r,sigmas(:,3,3),'Color',[0.6,0,1]),hold on;

    legend('11','12','13','21','22','23','31','32','33');
    xlabel('r[m]');
    ylabel('sigma');
    title('Conductivity tensor as function of k_x');
    drawnow;


    figure();
    plot(r,sigmas_dx(:,1,1),'Color',[1,0  ,0]),hold on;
    plot(r,sigmas_dx(:,1,2),'Color',[1,0.3,0]),hold on;
    plot(r,sigmas_dx(:,1,3),'Color',[1,0.6,0]),hold on;

    plot(r,sigmas_dx(:,2,1),'Color',[0.6,0.8,0  ]),hold on;
    plot(r,sigmas_dx(:,2,2),'Color',[0  ,0.8,0  ]),hold on;
    plot(r,sigmas_dx(:,2,3),'Color',[0  ,0.8,0.6]),hold on;

    plot(r,sigmas_dx(:,3,1),'Color',[0  ,0,1]),hold on;
    plot(r,sigmas_dx(:,3,2),'Color',[0.3,0,1]),hold on;
    plot(r,sigmas_dx(:,3,3),'Color',[0.6,0,1]),hold on;

    title('First derivative w.r.t. k_x of the conductivity tensor');
    xlabel('r[m]');
    ylabel('d/dk_x sigma');
    legend('11','12','13','21','22','23','31','32','33');
    drawnow;


    figure();
    plot(r,sigmas_dxx(:,1,1),'Color',[1,0  ,0]),hold on;
    plot(r,sigmas_dxx(:,1,2),'Color',[1,0.3,0]),hold on;
    plot(r,sigmas_dxx(:,1,3),'Color',[1,0.6,0]),hold on;

    plot(r,sigmas_dxx(:,2,1),'Color',[0.6,0.8,0  ]),hold on;
    plot(r,sigmas_dxx(:,2,2),'Color',[0  ,0.8,0  ]),hold on;
    plot(r,sigmas_dxx(:,2,3),'Color',[0  ,0.8,0.6]),hold on;

    plot(r,sigmas_dxx(:,3,1),'Color',[0  ,0,1]),hold on;
    plot(r,sigmas_dxx(:,3,2),'Color',[0.3,0,1]),hold on;
    plot(r,sigmas_dxx(:,3,3),'Color',[0.6,0,1]),hold on;
    title('Second derivative w.r.t. k_x of the conductivity tensor');
    xlabel('r[m]');
    ylabel('d^2/dk_x^2 sigma');
    legend('11','12','13','21','22','23','31','32','33');
    drawnow;






end

