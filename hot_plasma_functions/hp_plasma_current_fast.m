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

function [J_p] = hp_plasma_current_fast(r,E,k,curve_k_parr,omega,B,m,n,q,T,Jant)
%HP_PLASMA_CURRENT Summary of this function goes here
%   Calculates the plasma current in a fast fashion

    constants;
    %%%% 
    N_orders      = 4 ;
    N_per_Nk      = 20; %number of gridpoints per k-space gridpoints. Only has an effect when spine_interp = 1;
    spline_interp = 1 ; %use spline interpolation.
    %%%%
    
    N       = length(r);
    Nk_red  = 50;
    x_R     = max(r);
    x_L     = min(r);
    r0      = (x_R+x_L)/2;
    L       = x_R-x_L;
    x       = r-x_L;
    h       = L/(N-1);
    L       = L+h;
    
    if(~spline_interp)
        Nk_red = N;        
    end
    
    kx      = 1*ones(N,1);
    ky      = k(:,2);
    kz      = k(:,3); 
    
    Nk                      = floor(N/2);    
    kx_compressed_lowres    = linspace(-1,1,Nk_red)';
    kx_compressed_highres   = linspace(-1,1,N)';
    kxmax                   = 2*pi*Nk/L;    
    kx_lowres               = kx_compressed_lowres*kxmax;
    kx_highres              = kx_compressed_highres*kxmax;
    
    epsilons = zeros(Nk_red,9*N);
        
    %load Nk points of the dielectric tensor
    for ind = 1:Nk_red
        kx_test             = kx_lowres(ind)*ones(N,1);
        k_test              = [kx_test ky kz];
        epsilon             = hp_dielectric_tensor_swanson_per_species(N_orders,omega,k_test,B,T,q,n,m); 
        base = 0;
        for dim1 = 1:3
            for dim2 = 1:3
                epsilons(ind,base+(1:N))=epsilon(:,dim1,dim2);
                base = base+N;
            end
        end
    end

    %interpolate the dielectric tensor    
    epsilon_interp = interp1(kx_lowres,epsilons,kx_highres,'spline');
    
    %electric field spectrum
    E_spectrum     = fftshift(fft(E),1)/N;
    
    %allocate matrices for J_spectrum. 
    [J_spectrum_x,J_spectrum_y,J_spectrum_z] = deal(zeros(N,N));
    
    %calculate the (localized) spectrum of the plasma current
    for dim = 1:3
        J_spectrum_x = J_spectrum_x + epsilon_interp(:,(0+dim-1)*N+(1:N)).*E_spectrum(:,dim);
        J_spectrum_y = J_spectrum_y + epsilon_interp(:,(3+dim-1)*N+(1:N)).*E_spectrum(:,dim);
        J_spectrum_z = J_spectrum_z + epsilon_interp(:,(6+dim-1)*N+(1:N)).*E_spectrum(:,dim);
    end
    
    %Calculate the local plasma current using an IFFT
    Jp_local_x = N*diag(ifft(J_spectrum_x)).*exp(1i*kx_highres(1)*x);
    Jp_local_y = N*diag(ifft(J_spectrum_y)).*exp(1i*kx_highres(1)*x);
    Jp_local_z = N*diag(ifft(J_spectrum_z)).*exp(1i*kx_highres(1)*x);
    
    J_p = [Jp_local_x Jp_local_y Jp_local_z];
    
    J_p = J_p/(1i/(omega*eps0));
 
end

