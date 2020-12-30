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

function [sigma_stix] = hp_dielectric_tensor_stix(N,omega,k,B,T,q,n,m)
%SIGMA_STIX Summary of this function goes here
%   Detailed explanation goes here

constants;
%N        = 10;
e_x      = [1 0 0];
%calculate coordinate system stuff


%normalized magnetic field vector
norm_B   = sqrt(sum(B.^2,2));
e_b      = B./norm_B;

%stix frame
e_alpha  = e_x - sum(e_x.*e_b).*e_b;
e_alpha  = e_alpha./sqrt(sum(e_alpha.^2,2));
e_beta   = cross(e_b,e_alpha);

%wave vector coefficients
k_par    = sum(k.*e_b    ,2);
k_alpha  = sum(k.*e_alpha,2);
k_beta   = sum(k.*e_beta ,2);
k_perp   = sqrt(k_alpha.^2+k_beta.^2);

%plasma properties, as named in AORSA
alpha_th = thermal_velocity(T,m);
Omega    = cyclotron_frequency(q,norm_B,m);
rho      = larmour_radius(T,m,q,norm_B);
omega_p  = plasma_frequency(n,q,m);


Gamma    = 1/2*(k_perp.*rho).^2;
expGamma = exp(-Gamma);

[bessels, d_bessels]= load_all_bessels(N,Gamma);

sigma    = zeros(length(B),6);

kappa    = abs(k_par).*alpha_th;    %a sort-of parallel angular frequency


for l = -N:N

    zeta        = (omega-l*Omega)./kappa;
    
    %use approximation for the Dawson function
    %Zl          = 2/sqrt(pi)*dawson(zeta);  
    %takes the most time atm
    %Zl          = 2/sqrt(pi)*plasma_dispersion_function(zeta);
    
    %factor of 4 improvement can easily be made by "mirroring" and 
    %precomputing the value for every other l
    %[Il,dIl]    = besseli2(l,Gamma);
    Il =    bessels(:,abs(l)+1);
    dIl = d_bessels(:,abs(l)+1);
    
    
    if(sum(abs(k_par))==0)
        Pl          = omega_p.^2./(omega-l*Omega).*-1;       
        d_Pl        = zeros(size(kappa));
        if(sum(abs(k_perp))==0)
            sigma(:,1)  = sigma(:,1) +       0  ;
            if(abs(l)==1)
                sigma(:,2)  = sigma(:,2) +  1/2*Pl  ;
                sigma(:,3)  = sigma(:,3) + -1/2*Pl  ;
            end
            sigma(:,4)  = sigma(:,4) +       0  ;
            sigma(:,5)  = sigma(:,5) +       0  ;
            sigma(:,6)  = sigma(:,6) +       0  ; 
        else
            sigma(:,1)  = sigma(:,1) +       (Il-dIl).*Pl  ;
            sigma(:,2)  = sigma(:,2) +l^2  .*(Il    ).*Pl  ./Gamma;
            sigma(:,3)  = sigma(:,3) +l    .*(Il-dIl).*Pl  ;
            sigma(:,4)  = sigma(:,4) +                  0  ;
            sigma(:,5)  = sigma(:,5) +l    .*(Il    ).*d_Pl./Gamma;
            sigma(:,6)  = sigma(:,6) +       (Il-dIl).*d_Pl; 
        end
        
    else
        Zl          = plasma_dispersion_function(zeta);
        d_Zl        = -2*(1+zeta.*Zl); 

        Pl          = omega_p.^2./kappa.*Zl;       
        d_Pl        = omega_p.^2./kappa.*d_Zl;
        
        sigma(:,1)  = sigma(:,1) +       (Il-dIl).*Pl  ;
        sigma(:,2)  = sigma(:,2) +l^2  .*(Il    ).*Pl  ./Gamma;
        sigma(:,3)  = sigma(:,3) +l    .*(Il-dIl).*Pl  ;
        sigma(:,4)  = sigma(:,4) +zeta .*(Il    ).*d_Pl;
        sigma(:,5)  = sigma(:,5) +l    .*(Il    ).*d_Pl./Gamma;
        sigma(:,6)  = sigma(:,6) +       (Il-dIl).*d_Pl;        
    end
    
    
    
%     sigma(:,1)  = sigma(:,1) -1i  .*eps0*rho^2.*     (Il-dIl)       .*Pl  .*expGamma;
%     sigma(:,2)  = sigma(:,2) -1i  .*eps0      .*l^2.*(Il    )./Gamma.*Pl  .*expGamma;
%     sigma(:,3)  = sigma(:,3) -1   .*eps0      .*l  .*(Il-dIl)       .*Pl  .*expGamma;
%     sigma(:,4)  = sigma(:,4) +1i  .*eps0*zeta .*     (Il    )       .*d_Pl.*expGamma;
%     sigma(:,5)  = sigma(:,5) +1i/2.*eps0*rho  .*l  .*(Il    )./Gamma.*d_Pl.*expGamma;
%     sigma(:,6)  = sigma(:,6) + 1/2.*eps0*rho  .*     (Il-dIl)       .*d_Pl.*expGamma;



end
sigma = sigma.*eps0.*expGamma;
sigma(:,1)  = sigma(:,1).* -1i   .*rho.^2       ;
sigma(:,2)  = sigma(:,2).* -1i           ;
sigma(:,3)  = sigma(:,3).* -1                   ;
sigma(:,4)  = sigma(:,4).* +1i                  ;
sigma(:,5)  = sigma(:,5).* +1i/2 .*rho   ;
sigma(:,6)  = sigma(:,6).* + 1/2 .*rho          ;


sigma_stix = zeros(length(B),3,3);
sigma_stix(:,1,1) =  sigma(:,2)          + sigma(:,1).*k_beta.^2;
sigma_stix(:,1,2) =  sigma(:,3)          - sigma(:,1).*k_beta.*k_alpha;
sigma_stix(:,1,3) =  sigma(:,5).*k_alpha + sigma(:,6).*k_beta;

sigma_stix(:,2,1) = -sigma(:,3)          - sigma(:,1).*k_beta.*k_alpha;
sigma_stix(:,2,2) =  sigma(:,2)          + sigma(:,1).*k_alpha.^2;
sigma_stix(:,2,3) =  sigma(:,5).*k_beta  - sigma(:,6).*k_alpha;

sigma_stix(:,3,1) =  sigma(:,5).*k_alpha - sigma(:,6).*k_beta;
sigma_stix(:,3,2) =  sigma(:,5).*k_beta  + sigma(:,6).*k_alpha;
sigma_stix(:,3,3) =  sigma(:,4);

sigma_stix = sigma_stix*1i;

end

