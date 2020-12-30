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

function [y] = parabolic_profile(x,ap,y0,y_edge,alpha,lambda)
%PARABOLIC_PROFILE calculate parabolic profile with exponential decay in
% the edge, parabolic behaviour as copied from TOMCAT-U paper,
% exponential decay copied from cold plasma simulator
%   x: from r0-a to r0+a, a is the machine minor radius, r0 the major
%   radius
%   ap: plasma edge location (ap<a)
%   y0: density/temperature in the center
%   y_edge: density/temperature at the edge
%   alpha: steepness of parabolic profile
%   lambda: steepness of decay in the edge

    %gridpoints that are in the edge
    in_edge = abs(x)>ap;
    
    %central behaviour: parabola
    y=(y0-y_edge)*(1-(x/ap).^2).^alpha + y_edge;  
    
    %derivative at right edge (evaluate slightly left of ap)
    x_edge  = ap-1e-6;
    dy_edge =(y0-y_edge)*alpha*(1-(x_edge/ap)^2).^(alpha-1)*(-2*x_edge/ap^2);
    
    
    
    %y_in_edge   =    y_edge*exp(-b*(x-ap)))
    %d_y_in_edge = -b*y_edge*exp(-b*(x-ap)))
    
    b = -(dy_edge/y_edge);
    
    
    
    %edge behaviour: exponential decay (according to model D.van Eester)
    y(in_edge)=y_edge*exp(-abs(abs(x(in_edge))-ap)/lambda);
    
    %edge behavior:  exponential decay with smooth transition
    %y(in_edge) = y_edge*exp(-b*(abs(x(in_edge))-ap));
% %     
%     sigma = ap/3;
%     mu    = 0;
% %     
%      y = y0*exp(-((x-mu).^2)/(2*sigma)^2).^alpha + y_edge;
    
end

