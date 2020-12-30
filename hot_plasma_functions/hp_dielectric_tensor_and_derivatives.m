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
function [epsilon,d_epsilon,H_epsilon] = hp_dielectric_tensor_and_derivatives(N_orders,omega,k_test,dk,B,T,q,n,m)
%HP_DIELECTRIC_TENSOR_AND_DERIVATIVES Summary of this function goes here
%   Detailed explanation goes here
    N         = max(size(k_test,1),size(n,1));
    %steps for numerical differentiation in x,y,z. Step size is assumed to
    %be equal for all directions
    dkx         = [1 0 0]*dk;
    dky         = [0 1 0]*dk;
    dkz         = [0 0 1]*dk;
    
    %the
    epsilon     = hp_dielectric_tensor_swanson(N_orders,omega,k_test         ,B,T,q,n,m);
    

    %probe in x,y,z only
    epsilon_xp  = hp_dielectric_tensor_swanson(N_orders,omega,k_test+dkx     ,B,T,q,n,m);
    epsilon_yp  = hp_dielectric_tensor_swanson(N_orders,omega,k_test+dky     ,B,T,q,n,m);   
    epsilon_zp  = hp_dielectric_tensor_swanson(N_orders,omega,k_test+dkz     ,B,T,q,n,m);

    %probe twice in x,y,z only
    epsilon_xpp = hp_dielectric_tensor_swanson(N_orders,omega,k_test+2*dkx   ,B,T,q,n,m);
    epsilon_ypp = hp_dielectric_tensor_swanson(N_orders,omega,k_test+2*dky   ,B,T,q,n,m);    
    epsilon_zpp = hp_dielectric_tensor_swanson(N_orders,omega,k_test+2*dkz   ,B,T,q,n,m);
   
    %probe in 2 directions
    epsilon_xpyp = hp_dielectric_tensor_swanson(N_orders,omega,k_test+dkx+dky,B,T,q,n,m);
    epsilon_xpzp = hp_dielectric_tensor_swanson(N_orders,omega,k_test+dkx+dkz,B,T,q,n,m);
    epsilon_ypzp = hp_dielectric_tensor_swanson(N_orders,omega,k_test+dky+dkz,B,T,q,n,m);

    %derivatives with respect to x, y and z only
    epsilon_dx      = (-epsilon+epsilon_xp)/dk;
    epsilon_dy      = (-epsilon+epsilon_yp)/dk;
    epsilon_dz      = (-epsilon+epsilon_zp)/dk;
    
%     epsilon_dx      = (-3/2*epsilon+2*epsilon_xp-1/2*epsilon_xpp)/dk;
%     epsilon_dy      = (-3/2*epsilon+2*epsilon_yp-1/2*epsilon_ypp)/dk;
%     epsilon_dz      = (-3/2*epsilon+2*epsilon_zp-1/2*epsilon_zpp)/dk;
%     
    
    %double derivatives with respect to x, y and z only
    epsilon_dxx     = (epsilon-2*epsilon_xp+epsilon_xpp)/(dk^2);
    epsilon_dyy     = (epsilon-2*epsilon_yp+epsilon_ypp)/(dk^2);
    epsilon_dzz     = (epsilon-2*epsilon_zp+epsilon_zpp)/(dk^2);
    
    %mixed derivatives
    epsilon_dxy     = (-epsilon + epsilon_xpyp)/(dk*sqrt(2));
    epsilon_dxz     = (-epsilon + epsilon_xpzp)/(dk*sqrt(2));
    epsilon_dyz     = (-epsilon + epsilon_ypzp)/(dk*sqrt(2));
    
    %mixed derivatives (symmetric)
    epsilon_dyx     = epsilon_dxy;    
    epsilon_dzx     = epsilon_dxz;
    epsilon_dzy     = epsilon_dyz;
    
    %group the gradients of sigmas in one variable
    %first dimension: Nx1 (contribution per gridpoint)
    %second and third dimension: 3x3 (contribution per tensor element: epsxx .. epszz)
    %fourth dimension: 3x1 (contribution per derivative : d/dkx .. d/dkz)
    d_epsilon            = zeros(N,3,3,3);
    d_epsilon(:,:,:,1)   = epsilon_dx;
    d_epsilon(:,:,:,2)   = epsilon_dy;
    d_epsilon(:,:,:,3)   = epsilon_dz;
    
    %group the Hessian of sigma in one variable
    %first dimension: Nx1 (contribution per gridpoint)
    %second and third dimension: 3x3 (contribution per tensor element: epsxx .. epszz)
    %fourth and fifth dimension: 3x3 (contribution per derivative : d^2/dkx^2 .. d^2/dkz^2)
    H_epsilon            = zeros(N,3,3,3,3);
    H_epsilon(:,:,:,1,1) = epsilon_dxx;
    H_epsilon(:,:,:,1,2) = epsilon_dxy;
    H_epsilon(:,:,:,1,3) = epsilon_dxz;
    
    H_epsilon(:,:,:,2,1) = epsilon_dyx;
    H_epsilon(:,:,:,2,2) = epsilon_dyy;
    H_epsilon(:,:,:,2,3) = epsilon_dyz;
    
    H_epsilon(:,:,:,3,1) = epsilon_dzx;
    H_epsilon(:,:,:,3,2) = epsilon_dzy;
    H_epsilon(:,:,:,3,3) = epsilon_dzz;
    
    clearvars -except epsilon d_epsilon H_epsilon;
end

