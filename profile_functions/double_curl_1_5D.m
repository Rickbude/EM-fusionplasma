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

function [D0,D1,D2] = double_curl_1_5D(N,r,k,curve_k_parr)
%DOUBLE_CURL_1_5D Calculates the curl under the assumption that 
%E(r) = E(x)*exp(ikx*x+iky*y) The result is split in submatrices
%for constant, d/dx and d^2/dx^2 components.

    %wave number in y and z
    ky              = k(:,2);
    kz              = k(:,3); 

    D0              = zeros(N,3,3); %constant
    D1              = zeros(N,3,3); %contributions for 1st derivatives d/dx
    D2              = zeros(N,3,3); %contributions for 2nd derivatives d^2/dx^2
    
    %contribitions without derivatives with respect to x    
    D0(:,1,1)       = -(ky.^2 + kz.^2);
    D0(:,2,2)       = -kz.^2;
    D0(:,3,3)       = -ky.^2;
    D0(:,2,3)       = ky.*kz;
    D0(:,3,2)       = ky.*kz;
    
    %contributions with first derivatives with respect to x
    D1(:,1,2)       = -1i*ky;
    D1(:,2,1)       = -1i*ky;
    D1(:,1,3)       = -1i*kz;  
    D1(:,3,1)       = -1i*kz;
    
    %contributions with second derivatives with respect to x
    D2(:,2,2)       = 1;
    D2(:,3,3)       = 1;
    
    %extra terms due to curvature
    if(curve_k_parr)
        D0(:,3,1)      = D0(:,3,1) + 1i*kz./r; %should be +1i*kz./r
        D0(:,1,3)      = D0(:,1,3) - 1i*kz./r; %should be -1i*kz./r
        D0(:,3,3)      = D0(:,3,3) - 1./(r.^2);
        D0(:,2,1)      = D0(:,2,1) - 1i*ky./r;
        
        D1(:,2,2)      = D1(:,2,2) + 1./r;
        D1(:,3,3)      = D1(:,3,3) + 1./r;
    end
    
    clearvars -except D0 D1 D2;
end

