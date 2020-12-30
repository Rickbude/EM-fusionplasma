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

function [Hessian_r0,Hessian_r1,Hessian_r2] = Hessian_r_1_5D(N,k)
%HESSIAN_R_1_5D Summary of this function goes here
%   Detailed explanation goes here

    %wave number in y and z
    ky              = k(:,2);
    kz              = k(:,3); 

    Hessian_r0 = zeros(N,3,3);
    Hessian_r1 = zeros(N,3,3);
    Hessian_r2 = zeros(N,3,3);
    
    Hessian_r0(:,2,2) = -ky.^2 ;
    Hessian_r0(:,2,3) = -ky.*kz;
    Hessian_r0(:,3,2) = -ky.*kz;
    Hessian_r0(:,3,3) = -kz.^2 ; 
    
    Hessian_r1(:,1,2) = 1i*ky;
    Hessian_r1(:,2,1) = 1i*ky;
    Hessian_r1(:,1,3) = 1i*kz;
    Hessian_r1(:,3,1) = 1i*kz; 
    
    Hessian_r2(:,1,1) = 1;
    
    clearvars -except Hessian_r0 Hessian_r1 Hessian_r2;
end

