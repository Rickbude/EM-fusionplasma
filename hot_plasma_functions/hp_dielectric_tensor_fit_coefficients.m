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

function [ak] = hp_dielectric_tensor_fit_coefficients(profiles,flags,settings)
%HP_DIELECTRIC_TENSOR_FIT_COEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here
%CP_CALC_ELECTRIC_FIELD 
%   Detailed explanation goes here
    constants;
    
    %extract r and r-dependent profiles
    B            = profiles.Bb;
    k            = profiles.k;
    m            = profiles.m;
    n            = profiles.n;
    q            = profiles.q;
    T            = profiles.T;
    r            = profiles.r;
    %extract simulation settings
    omega        = settings.omega;
    N            = settings.N;
    N_orders     = settings.N_harmmax;          %number of terms in the dielectric tensor sum   
    Nk           = settings.Nk;                 %number of k-space sample points
    Np           = settings.Np;                 %polynomial order
    zeta         = settings.zeta;               %k_perp*rho_L number for throughout domain
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    kx     = 1*ones(N,1);
    ky     = k(:,2);
    kz     = k(:,3); 
    
    fprintf('Start calculating dielectric tensor (%.0f ms)\n',toc*1000);

    %calculate the larmour radius of each species
    rho_L     = larmour_radius(T,m,q,B);
    %find the smallest larmour radius, excluding the electron radius
    rho_L_min = min(rho_L(:,2:end),[],2);
    
    kxmax = zeta./rho_L_min;
    
    kx_compressed   = linspace(-1,1,Nk)';
    kx              = kx_compressed;
    
    ak = zeros(N,3,3,Np+1); 
    
    w      = waitbar(0,'loading tensor');
    tprev  = toc;
    
    epsilons = zeros(N,Nk,3,3);
    for ind = 1:Nk
        kx_test             = kx_compressed(ind)*kxmax;
        k_test              = [kx_test ky kz];
        epsilons(:,ind,:,:) = hp_dielectric_tensor_swanson(N_orders,omega,k_test,B,T,q,n,m); 
        if(toc-tprev > 0.5)
            waitbar(ind/N,w,sprintf('loading tensor (%i / %i)',ind,Nk));
            tprev = toc;
        end
    end
    close(w)
    fprintf('Start fitting polynomials through tensor components(%.0f ms)\n',toc*1000);
    w      = waitbar(0,'calculating coefficients');
    errors  = zeros(N,3,3);
    
    %set which epsilon components (eps_xx, eps_xy, etc) take also negative
    %kx into account. This is most suitable for the smooth behaving
    %components around kx = 0. 110 110 001 is also a potential setting
    smooth = [
        1 1 1
        1 1 1
        1 1 1       
    ];

    %construct the Vandermonde matrix only once, used for finding the
    %coefficients. This is an advantage over repeated calles to polyfit,
    %where every iteration this (same) matrix needs to be reconstructed
    %MATLAB's built-in polyfit achieves the same result. it is much slower,
    %but provides warnings when the fit is badly conditioned. When 
    %examining poor numerical results, consider enabling polyfit again.
    kx_compressed_pos = kx_compressed(kx>0);
    kx_matrix         = zeros(length(kx),Np+1);
    kx_matrix_pos     = zeros(nnz(kx>0) ,Np+1);
    
    %construct the Vandermonde matrix
    for i = 1:Np+1
        kx_matrix(:,i)    = kx_compressed.^(Np-i+1);
        kx_matrix_pos(:,i)= kx_compressed_pos.^(Np-i+1);
    end
    
    %a decent amount of time is saved by calculating this (pseudo-)inverse
    %only once, and using it later to calculate the fit
    kx_matrix_inv     = pinv(kx_matrix);
    kx_matrix_pos_inv = pinv(kx_matrix_pos);
    
    %"decompression coefficients"
    decompression = kxmax.^fliplr(0:Np).';
            
    do_plot = 0;
    counter = 1;
    for dim1 = 1:3
        for dim2 = 1:3
            epsilon = squeeze(epsilons(:,:,dim1,dim2)).';
            for ind = 1:N
                
                compression = kxmax(ind);
                kx = kx_compressed*compression;
        
                epsilon_dimdim = epsilon(:,ind);                
                
                if(nnz(kx_compressed<0)==0 || smooth(dim1,dim2))                   
                   %p = polyfit(kx_compressed,epsilon_dimdim,Np).';
                   %can probably be sped up further using QR decomposition
                   %p_compressed = kx_matrix\epsilon_dimdim;   
                   %or, by calculating the inverse only once..
                    p_compressed  = kx_matrix_inv*epsilon_dimdim;
                else
                   %p = polyfit(kx_compressed(kx>0),epsilon_dimdim(kx>0),Np).';
                   %p_compressed = kx_matrix_pos\epsilon_dimdim(kx>0);
                    p_compressed  = kx_matrix_pos_inv*epsilon_dimdim(kx>0);
                end
                %reverse the kx-compression by compensating in the
                %coefficients
                
                p = p_compressed./decompression(:,ind);
                ak(ind,dim1,dim2,:) = p;      
                
                if(nnz(kx<0)==0 || smooth(dim1,dim2))      
                   %the curve can be reconstructed using polyval, but as
                   %the kx_matrix is available, the direct method is faster
                   %fitted = polyval(p,kx);
                   fitted = kx_matrix*p_compressed;
                   actual = epsilon_dimdim;
                   %errors(ind,dim1,dim2) = sum(abs((fitted-actual)./(actual)))/Nk;        
                   errors(ind,dim1,dim2) = sum(abs((fitted-actual)))/Nk;      
                else
                    %fitted = polyval(p,kx(kx>0));
                    fitted = kx_matrix_pos*p_compressed;
                    actual = epsilon_dimdim(kx>0);
                    errors(ind,dim1,dim2) = sum(abs((fitted-actual)./(actual)))/Nk;
                end      
                
               
                if(mod(ind,floor(N/20))==0 && do_plot)
                    subplot(3,3,counter);
                    plot(kx,real(epsilon_dimdim)),hold on;
                    plot(kx,real(polyval(p,kx)),'-.'),hold off;
                    drawnow;
                    pause(0.1);
                end               
                
                if(toc-tprev > 0.5)
                    waitbar((counter+ind/N)/9,w,sprintf('calculating coefficients (%i / %i)',round((counter+ind/N)/9*N),N));
                    tprev = toc;
                end                
            end
            counter = counter + 1;
        end           
    end
    close(w);
    clearvars -except ak;
   
end

