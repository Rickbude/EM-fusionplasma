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

function [J_plasma] = hp_plasma_current_modified(E,profiles,flags,settings)
%HP_PLASMA_CURRENT_MODIFIED Summary of this function goes here
%   Detailed explanation goes here
    constants;
    
    %extract r and r-dependent profiles
    B            = profiles.Bb;
    k            = profiles.k;
    m            = profiles.m;
    n            = profiles.n;
    q            = profiles.q;
    T            = profiles.T;
    Jant         = profiles.Jant;
    r            = profiles.r;
    %extract simulation settings
    omega        = settings.omega;
    N            = settings.N;
    N_orders     = settings.N_harmmax;          %number of terms in the dielectric tensor sum   
    x_w          = settings.x_wall;             %wall location
    Nk           = settings.Nk;                 %number of k-space sample points
    Np           = settings.Np;                 %polynomial order
    zeta         = settings.zeta;               %k_perp*rho_L number for throughout domain
    %extract simulation flags
    curve_k_parr = flags.curve_k_parr;          %curved (1) or straight (0) magnetic field lines 
    aperature    = flags.aperature;             %use aperature-type excitation (1) or current sheet (0)
    
    [~,coeff,~] = fd_coefficients(Np);
    ncoeff     = length(coeff);
    s           = -floor(ncoeff/2):floor(ncoeff/2);  
    E_derivative = zeros(N,3,Np+1);
    for p = 0:Np
        c = coeff(p+1,:);
        for i=1:ncoeff
            E_derivative(:,:,p+1) = E_derivative(:,:,p+1) + c(i)*circshift(E,s(i),1);
        end
    end
    
    ak = hp_dielectric_tensor_fit_coefficients(profiles,flags,settings);
    ak(:,1,1,1) = ak(:,1,1,1) - 1;
    ak(:,2,2,1) = ak(:,2,2,1) - 1;
    ak(:,3,3,1) = ak(:,3,3,1) - 1;
    
    
    J_plasma_3x3 = zeros(N,3,3);
    for dim1 = 1:3
        for dim2 = 1:3
            for p = 1:Np+1
                J_plasma_3x3(:,dim1,dim2) = J_plasma_3x3(:,dim1,dim2) + (1i)^(-(p-1))*squeeze(ak(:,dim1,dim2,end-p+1)).*E_derivative(:,dim2,p);                
            end
        end
    end
    
    
    figure(89)
    counter = 1;
    for i = 1:Np+1
       subplot(Np+1,3,counter+0);
       plot(r,real(E_derivative(:,1,i))),hold on;
       plot(r,imag(E_derivative(:,1,i))),hold on;
       subplot(Np+1,3,counter+1);
       plot(r,real(E_derivative(:,2,i))),hold on;
       plot(r,imag(E_derivative(:,2,i))),hold on;
       subplot(Np+1,3,counter+2);
       plot(r,real(E_derivative(:,3,i))),hold on;
       plot(r,imag(E_derivative(:,3,i))),hold on;
       counter= counter+3;
    end
    
    figure(90)
    counter = 1;
    J_plasma = zeros(N,3);
    for dim1 = 1:3
        for dim2 = 1:3
            subplot(3,3,counter);
            plot(r,real(J_plasma_3x3(:,dim1,dim2)));
            J_plasma(:,dim1) = J_plasma(:,dim1)+J_plasma_3x3(:,dim1,dim2);
            counter = counter+1;
        end
        
    end
    J_plasma = J_plasma/(1i/(omega*eps0));
    
end

