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

function [sigma_swanson] = hp_dielectric_tensor_swanson_per_species(N_orders,omega,k,B,Tj,qj,nj,mj)
%SIGMA_STIX Summary of this function goes here
%   Detailed explanation goes here

%assumption:
%k = kx ky kz
%B = |B|


N        = max(size(nj,1),size(k,1));

kx       = k(:,1);
ky       = k(:,2);
kz       = k(:,3);

norm_B   = sqrt(sum(B.^2,2));

k_perp   = sqrt(kx.^2 + ky.^2);
k_parr   = kz;

%plasma frequency, cyclotron frequency, Larmor radius and
%thermal velocity
omega_p  = plasma_frequency(nj,qj,mj);
Omega_c  = cyclotron_frequency(qj,norm_B,mj);
rho_L    = larmour_radius(Tj,mj,qj,norm_B);
v_th     = thermal_velocity(Tj,mj);

%dimensionless
lambda   = 1/2*(k_perp.*rho_L).^2;


[K0,K1,K2,K3,K4,K5] = deal(zeros(N,1));

k_parr_limit = sum(abs(k_parr))==0;
k_perp_limit = sum(abs(k_perp))==0;


if(k_perp_limit)
    bessels_exp                  = zeros(N,N_orders+1);
    d_bessels_exp                = zeros(N,N_orders+1);
    
    bessels_exp  (:,1)  = 1;    %only 1   at n =  0
    d_bessels_exp(:,2)  = 1/2;  %only 1/2 at n = +- 1
else    
    [bessels_exp, d_bessels_exp] = load_all_bessels(N_orders,lambda);
end


% if(k_parr_limit && k_perp_limit)
%     %TODO: both k_parr and k_perp are 0
%     disp('Both k_perp = 0 and k_parr = 0');
%     %return;
% elseif(k_parr_limit)
%     disp('Only k_parr = 0');
% elseif(k_perp_limit)
%     %TODO: only k_perp = zero
%     disp('Only k_perp = 0 ');
%     %return;        
% else
%     %disp('Both k_perp ~= 0 and k_parr ~= 0 ');
% end


for n = -N_orders:N_orders
    In  =    bessels_exp(:,abs(n)+1);
    dIn =  d_bessels_exp(:,abs(n)+1);
    
    if(k_parr_limit && k_perp_limit)
        %TODO: both k_parr and k_perp are 0
         Zn = -1;
        dZn =  1; 
        
        %prod1       = omega_p.^2.*exp(-lambda)./(omega);
        prod1       = omega_p.^2./(omega);
        prod2       = v_th./(omega+n*Omega_c);
        
        fact1       = prod1.*prod2./v_th;
        fact2       = prod1.*prod2.*k_perp./(2*Omega_c);

        K0  = K0 + 2* fact1.*lambda.*(In-dIn).*Zn;
        %handle limit of besseli/lambda
        if(n == 1 || n == -1)
           %K1  = K1 +    fact1.*l^2.*Il .*Zl./lambda;
            K1  = K1 +    fact1.*n^2.*dIn.*Zn;
        end  
        
        K2  = K2 + 1i*fact1.*sign(qj).*n.*(In-dIn).*Zn;
        %K3  = K3 + 0;
        K3  = K3 + -1*fact1.*In;
        %handle limit of besseli/lambda
        if(n == 1 || n == -1)
           %K4  = K4 +    fact2.*l.*Il .*dZl./lambda;
           % K4  = K4 +    fact2.*l.*dIl.*dZl;
           K4 = K4+0;
        end
        K5  = K5 + 1i*fact2.*sign(qj).*(In-dIn).*dZn;
    elseif(k_parr_limit)
        %TODO: only k_parr = 0
         Zn = -1;
        dZn =  1; 
        
        %prod1       = omega_p.^2.*exp(-lambda)./(omega);
        prod1       = omega_p.^2./(omega);
        prod2       = v_th./(omega+n*Omega_c);
        
        fact1       = prod1.*prod2./v_th;
        fact2       = prod1.*prod2.*k_perp./(2*Omega_c);

        K0  = K0 + 2* fact1.*lambda.*(In-dIn).*Zn;
        K1  = K1 +    fact1.*n^2.*In .*Zn./lambda;
        K2  = K2 + 1i*fact1.*sign(qj).*n.*(In-dIn).*Zn;
        %K3  = K3 + 0;
        K3  = K3 + -1*fact1.*In;
        K4  = K4 + 0;
        %K4  = K4 +    fact2.*l.*Il.*dZl;
        K5  = K5 + 1i*fact2.*sign(qj).*(In-dIn).*dZn;
    elseif(k_perp_limit)
        %both k_perp and k_parr are non-zero -> use full expressions
        %dimensionless
        zeta        = (omega./Omega_c+n)./(abs(k_parr).*rho_L);
        
        %plasma dispersion number and its derivative
        Zn  = plasma_dispersion_function(zeta);
        dZn = -2*(1+zeta.*Zn); 
        
       % prod1       = omega_p.^2.*exp(-lambda)./(omega.*kz);
        prod1       = omega_p.^2./(omega.*kz);
        
        fact1       = prod1./v_th;
        fact2       = prod1.*k_perp./(2*Omega_c);

        K0  = K0 + 2* fact1.*lambda.*(In-dIn).*Zn;
        %handle limit of besseli/lambda
        if(n == 1 || n == -1)
           %K1  = K1 +    fact1.*l^2.*Il .*Zl./lambda;
            K1  = K1 +    fact1.*n^2.*dIn.*Zn;
        end        
        K2  = K2 + 1i*fact1.*sign(qj).*n.*(In-dIn).*Zn;
        K3  = K3 + -1*fact1.*In.*zeta.*dZn;
        %handle limit of besseli/lambda
        if(n == 1 || n == -1)
           %K4  = K4 +    fact2.*l.*Il .*dZl./lambda;
            K4  = K4 +    fact2.*n.*dIn.*dZn;
        end
        K5  = K5 + 1i*fact2.*sign(qj).*(In-dIn).*dZn;        
    else
        %both k_perp and k_parr are non-zero -> use full expressions
        zeta        = (omega./Omega_c+n)./(abs(k_parr).*rho_L);
        %zeta         = (omega + n*Omega_c)./(kz.*v_th);

        Zn  = plasma_dispersion_function(zeta);
        dZn = -2*(1+zeta.*Zn); 
        
        %prod1       = omega_p.^2.*exp(-lambda)./(omega.*kz);
        prod1       = omega_p.^2./(omega.*kz);
        
        fact1       = prod1./v_th;
        fact2       = prod1.*k_perp./(2*Omega_c);

        K0  = K0 + 2* fact1          .*(In-dIn).*Zn .*lambda;
        K1  = K1 +    fact1          .*(In    ).*Zn .*n^2./lambda;
        K2  = K2 + 1i*fact1.*sign(qj).*(In-dIn).*Zn .*n;
        K3  = K3 + -1*fact1          .*(In    ).*dZn.*zeta;
        K4  = K4 +    fact2          .*(In    ).*dZn.*n./lambda;
        K5  = K5 + 1i*fact2.*sign(qj).*(In-dIn).*dZn;
        
        
    end
    
    

    
    
end



phi = atan2(ky,real(kx));
% phi = (1+sign(kx))/2*pi;
% phi = 0;

sigma_swanson = zeros(N,3,3);
sigma_swanson(:,1,1) =  K1 + sin(phi).^2.*K0;
sigma_swanson(:,1,2) =  K2 - cos(phi).*sin(phi).*K0;
sigma_swanson(:,1,3) =  cos(phi).*K4 + sin(phi).*K5;

sigma_swanson(:,2,1) = -K2 - cos(phi).*sin(phi).*K0;
sigma_swanson(:,2,2) =  K1 + K0.*cos(phi).^2;
sigma_swanson(:,2,3) = -K5.*cos(phi) + K4.*sin(phi);

sigma_swanson(:,3,1) =  K4.*cos(phi) - K5.*sin(phi);
sigma_swanson(:,3,2) =  K5.*cos(phi) + K4.*sin(phi);
sigma_swanson(:,3,3) =  K3;


% sigma_swanson = zeros(N,3,3);
% sigma_swanson(:,1,1) =  K1;
% sigma_swanson(:,1,2) =  K2;
% sigma_swanson(:,1,3) =  K4;
% 
% sigma_swanson(:,2,1) = -K2;
% sigma_swanson(:,2,2) =  K0+K1;
% sigma_swanson(:,2,3) = -K5;
% 
% sigma_swanson(:,3,1) =  K4;
% sigma_swanson(:,3,2) =  K5;
% sigma_swanson(:,3,3) =  K3;

clearvars -except sigma_swanson;


end

