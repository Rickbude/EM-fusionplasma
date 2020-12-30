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

function [J_p] = hp_plasma_current(r,E,k,curve_k_parr,omega,B,m,n,q,T,Jant)
%HP_PLASMA_CURRENT Summary of this function goes here
%   Detailed explanation goes here

    %%%% 
    N_orders      = 5 ;
    N_per_Nk      = 20; %number of gridpoints per k-space gridpoints. Only has an effect when spine_interp = 1;
    spline_interp = 1 ; %use spline interpolation.
    %%%%
    
    N       = length(r);
    Nk_red  = max(ceil(N/N_per_Nk),20);
    x_R     = max(r);
    x_L     = min(r);
    r0      = (x_R+x_L)/2;
    L       = x_R-x_L;
    x       = r-r0;
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
    dk                      = (max(kx_highres)-min(kx_highres))/(N-1);
    
    
    
    
    epsilons = zeros(N,Nk_red,3,3);
    w      = waitbar(0,'loading tensor');
    tprev  = toc;
    for ind = 1:Nk_red
        kx_test             = kx_lowres(ind)*ones(N,1);
        k_test              = [kx_test ky kz];
        epsilons(:,ind,:,:) = hp_dielectric_tensor_swanson_per_species(N_orders,omega,k_test,B,T,q,n,m); 
        if(toc-tprev > 0.5)
            waitbar(ind/Nk_red,w,sprintf('Loading tensor(%i / %i)',ind,Nk_red));
            tprev = toc;
        end
    end
    close(w);

    
    w      = waitbar(0,'Interpolating');
    tprev  = toc;
    
    %interpolate the dielectric tensor
    epsilon_interp = zeros(N,N,3,3);
    counter = 1;    
    if(spline_interp)
        for dim1 = 1:3
            for dim2 = 1:3               
                %load dielectric tensor component{dim1,dim2}. contains 
                %values for all gridpoints and k-space sample points
                epsilon  = squeeze(epsilons(:,:,dim1,dim2));       
                %interp1 can handle interpolation for N vectors
                %simultaniously very efficiently
                e_interp = interp1(kx_lowres,epsilon(:,:).',kx_highres,'spline');
                epsilon_interp(:,:,dim1,dim2) = e_interp.';                
                if(toc-tprev > 0.5)
                    waitbar(counter/9,w,sprintf('interpolating (%i / %i)',counter,9));
                    tprev = toc;
                end                  
                counter = counter + 1;
            end            
        end
    else
        epsilon_interp = epsilons;
    end
    close(w)

    do_plot = 0;    
    for ind = 1:N
        counter = 1;
        for dim1 = 1:3
            for dim2 = 1:3
                if(do_plot && mod(ind,20)==0)
                    subplot(3,3,counter);
                    plot(kx_lowres ,abs(squeeze(epsilons(ind,:,dim1,dim2)))      ,'.'),hold on;
                    plot(kx_highres,abs(squeeze(epsilon_interp(ind,:,dim1,dim2))),'--'),hold on;
                    hold off;
                    drawnow;
                end
                counter = counter+1;
            end            
        end
    end    
    
    kx  = 2*pi/L*(-Nk:Nk)';  

    E_spectrum = fftshift(fft(E),1);

    J_p = zeros(N,3);
    x   = r-x_L;
    w      = waitbar(0,'Calculating plasma current');
    tprev  = toc;
    
%     J_spectrum = zeros(N,N,3);
%     for dim1 = 1:3
%         for dim2 = 1:3
%             epsilon = squeeze(epsilon_interp(:,:,dim1,dim2));
%             J_spectrum(:,:,dim1) = J_spectrum(:,:,dim1) + epsilon.*E_spectrum
%             for ind = 1:N
%                 J_spectrum(ind,:,dim1) = J_spectrum(ind,:,dim1) + 
%             end
%         end        
%     end
    
    for ind = 1:N
         fact = exp(1i*kx_highres*x(ind));
         epsilon  = squeeze(epsilon_interp(ind,:,:,:));         
         J_p_spectrum = zeros(N,3);
         for dim = 1:3
             J_p_spectrum(:,dim) = sum(E_spectrum.*squeeze(epsilon(:,dim,:)),2);
         end
         fact2     = exp(1i*kx_highres(1)*x);
         
         %Calculate the local plasma current using an IFFT
         %also works, but slightly slower
         %J_p_local = N*ifft(J_p_spectrum).*fact2;
         %J_p_local = J_p_local(ind,:);
         
         J_p_local  = sum(J_p_spectrum.*fact,1);         
         J_p(ind,:) = J_p_local;
         

         if(toc-tprev > 0.5)
            waitbar(ind/N,w,sprintf('Calculating plasma current(%i / %i)',ind,N));
            tprev = toc;
        end       
    end
    
    close(w);    
  
end

