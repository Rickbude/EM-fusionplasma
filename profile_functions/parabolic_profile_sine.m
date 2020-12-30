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
function [y] = parabolic_profile_sine(x,ap,y0,y_edge,alpha,lambda)
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
    
    w_edge = max(x)-ap; %edge width
    
    %central behaviour: parabola
    y=(y0-y_edge)*(1-(x/ap).^2).^alpha + y_edge;  
    
    %derivative at right edge (evaluate slightly left of ap)
    x_edge  = ap-1e-6;
    dy_edge =(y0-y_edge)*alpha*(1-(x_edge/ap)^2).^(alpha-1)*(-2*x_edge/ap^2);
    N       = length(x);
    dx      = (max(x)-min(x))/(N-1);
    y_in_edge = y(in_edge);
    dy_edge = (y_in_edge(2)-y_in_edge(1))/dx;
    
    
    y(in_edge) = 0;
    
    x_e = w_edge/2;    
    
    fun = @(b) dy_edge/(b*sin(b*x_e))*(1-cos(b*x_e))-y_edge;
    
    x0  = x_e;
    
    bmin = lsqnonlin(fun,x0);
    amin = y_edge/(1-cos(bmin*x_e));
    
    y(in_edge) = 0;
    
    in_sine = bitand(abs(x)>ap,abs(x)<ap+x_e);
    x_sine = x_e-(abs(x(in_sine))-ap);
    y(in_sine) = amin*(1-cos(bmin*x_sine));
    
    
    %%%%%%%%%%%% RESET %%%%%%%%%%%%%%%%%%%%
    ap_new = ap+w_edge/2;
    y  = y0*(1+cos(pi*x/ap_new));
    y(abs(x)>ap_new) = 0;
    y = y+y_edge;
    
end

