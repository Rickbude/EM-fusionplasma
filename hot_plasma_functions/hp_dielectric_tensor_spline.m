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

function [epsilons] = hp_dielectric_tensor_spline(profiles,flags,settings)
%HP_DIELECTRIC_TENSOR_SPLINE Summary of this function goes here
%   Detailed explanation goes here

    %load constants
    constants;
    
    %extract r and r-dependent profiles
    B            = profiles.Bb;
    k            = profiles.k;
    m            = profiles.m;
    n            = profiles.n;
    q            = profiles.q;
    T            = profiles.T;
    r            = profiles.r;
    
    %flags
    cold_plasma  = flags.cold_plasma;
    
    %extract simulation settings
    omega        = settings.omega;
    N            = settings.N;
    N_orders     = settings.N_harmmax;          %number of terms in the dielectric tensor sum   
    Nk_red       = settings.Nk;                 %number of k-space sample points
    
    %Don't sample the k-space more than needed
    if(Nk_red>N)
       fprintf('Only sampling k-space %i times instead of Nk times (%i)\n',N,Nk_red);
       Nk_red = N; 
    end

    x_R     = max(r);
    x_L     = min(r);
    L       = x_R-x_L;
    h       = L/(N-1);
    L       = L+h;
    
    ky      = k(:,2);
    kz      = k(:,3); 
    
    Nk                      = floor(N/2);    
    kx_compressed_lowres    = linspace(-1,1,Nk_red)';    
    kxmax                   = 2*pi*Nk/L;    
    kx_lowres               = kx_compressed_lowres*kxmax;
    
    
    epsilons = complex(zeros(Nk_red,N,3,3));
    
    w = waitbar(0,sprintf('loading dielectric tensor (%i /%i)',0,Nk_red));
    t_prev = toc;    
    %load Nk points of the dielectric tensor
    for ind = 1:Nk_red
        kx_test             = kx_lowres(ind)*ones(N,1);
        k_test              = [kx_test ky kz];
        if(cold_plasma)
            epsilon             = cp_dielectric_tensor(omega,B,m,n,q); 
        else
            epsilon             = hp_dielectric_tensor_swanson(N_orders,omega,k_test,B,T,q,n,m); 
        end
                
        epsilons(ind,:,:,:) = epsilon; 
        %limit the waitbar updates to 10 Hz. (too frequent updates give
        %significant slowdown)
        if(toc - t_prev > 0.1)
            waitbar(ind/Nk_red,w,sprintf('loading dielectric tensor (%i /%i)',ind,Nk_red));
            t_prev = toc;
        end
    end
    close(w);
    
    %interpolate intermediate k samples using spline interpolation. Skip
    %this step if it is not required.
    if(N ~= Nk_red)
        
        disp('Spline interpolating dielectric tensor');
        %reshape such that a single call to interp1 suffices
        epsilons = reshape(epsilons,[Nk_red N*9]);

        %interpolate the dielectric tensor    
        kx_compressed_highres   = linspace(-1,1,N)';
        kx_highres              = kx_compressed_highres*kxmax;    
        epsilon_interp          = interp1(kx_lowres,epsilons,kx_highres,'spline');

        %reshape back
        epsilons = reshape(epsilon_interp,[N N 3 3]);
    end
    
    clearvars -except epsilons;  
end

