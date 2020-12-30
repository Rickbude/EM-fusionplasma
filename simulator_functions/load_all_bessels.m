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
function [bessels,d_bessels] = load_all_bessels(N,Gamma)
%LOAD_ALL_BESSELS Summary of this function goes here
%   Detailed explanation goes here

    bessels  = zeros(length(Gamma),N+3);
    d_bessels= zeros(length(Gamma),N+3);
    
    %calculate every other besseli function
    for n = 0:2:N+2
        %bessels(:,n+1)   = besseli(n,Gamma);
        bessels(:,n+1)    = bessel_times_exponent_approx(n,Gamma);
    end

    %interpolate other bessels
    for n = 1:2:N+1
        I_min            = bessels(:,abs(n-1)+1);
        I_plus           = bessels(:,abs(n+1)+1);
        bessels(:,n+1)   = Gamma/(2*n).*(I_min-I_plus);   
    end

    %calculate derivatives 
    for n = 0:N+1
        I_min            = bessels(:,abs(n-1)+1);
        I_plus           = bessels(:,abs(n+1)+1);   
        d_bessels(:,n+1) = 1/2*(I_min+I_plus);
    end
    bessels   =  bessels(:,1:N+1);
    d_bessels =d_bessels(:,1:N+1);

end

