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

function [profiles,flags,settings] = load_profile_override(input_file,override)
%LOAD_PROFILE Summary of this function goes here
%   Detailed explanation goes here
    %run(strcat('plasma_profiles/',profile_name,'.m'));
    run(input_file);
    %override the loaded profile parameters by evaluating the override data
    keySet = fieldnames(override);
    Nkeys  = length(keySet);    
    for keyind = 1:Nkeys
       key = char(keySet(keyind));
       feval(@()assignin('caller',key,override.(key))); 
    end
    
    N
    
    
    constants; 
    
    %complex frequency
    omega   = 2*pi*f*(1+1i*nuom); 

    %derived machine parameters
    x_R     = r0+a;
    x_L     = r0-a;
    r_a     = r0+x_ant;
    
    %grid
    r       = linspace(x_L,x_R,N)';   

    Bbx     = 0 .*ones(size(r));
    Bby     = 0 .*ones(size(r));
    Bbz     = B0./(1+(r-r0)/r0);     %1/R dependence, based on magnetic field on axis B0

    Bb      = [Bbx Bby Bbz];     

    n = zeros(N,n_species+1);
    T = zeros(N,n_species+1);
    m = zeros(N,n_species+1);
    q = zeros(N,n_species+1);
    
    %%%% test continuous density - introduce new flag
    x = r-r0;
    
    if(iatanprofile)
        %iatan profile, copied from cold plasma simulator D. van Eester
        ne=n0*(atan(-(x-ap)/ap/fracap)/pi+0.5).*(atan((x+ap)/ap/fracap)/pi+0.5);              
    else
        ne = parabolic_profile(x,ap,n0,n0_edge,alpha_n,lambda_n);
    end
%     ne(x>(r_a-r0)) = 0;


    Te = parabolic_profile(x,ap,T0,T0_edge,alpha_T,lambda_T);
    Te = Te*e/kB;
    
    n(:,1) = ne;
    m(:,1) = me;
    q(:,1) = -e;
    T(:,1) = Te;
    
    if(~T_ion)
       T_ion = T0*ones(n_species); 
    end
    
    if(~alpha_n_ion)
       alpha_n_ion = alpha_n*ones(n_species);       
    end
    
    if(~alpha_T_ion)
       alpha_T_ion = alpha_T*ones(n_species);       
    end
        
    for i = 1:n_species   
        m(:,1+i)  = A_ion(i)*mp;           
        q(:,1+i)  = Z_ion(i)*e;
%         if(T_ion)           
        T(:,1+i)  = parabolic_profile(x,ap,T_ion(i),T0_edge,alpha_T_ion(i),lambda_T)*e/kB;
        n0_i      = conc_ion(i)*n0./Z_ion(i);
        n0_i_edge = conc_ion(i)*n0_edge./Z_ion(i);
        n(:,1+i)  = parabolic_profile(x,ap,n0_i,n0_i_edge,alpha_n_ion(i),lambda_n);
%         else
%             T(:,1+i) = Te;
%             n(:,1+i) = conc_ion(i)*ne./Z_ion(i);
%         end        
    end
    
    I_ant = [0 I_ant_y I_ant_z];
    Jant  = antenna_current(I_ant,r,r_a);

    if(curve_k_parr)
        disp('Using a variable k_z profile');
        k_parr  = ntor./r;
    else
        disp('Using a flat k_z profile');
        k_parr  = ones(N,1)*ntor/r0;
    end
    k       = [zeros(N,1) ky*ones(N,1)  k_parr];    
    
    %r and r-dependent profiles
    profiles.r          = r;
    profiles.n          = n;
    profiles.T          = T;
    profiles.m          = m;
    profiles.q          = q;
    profiles.Jant       = Jant;
    profiles.k          = k;
    profiles.Bb         = Bb;
    
    %flags
    flags.curve_k_parr  = curve_k_parr;
    flags.aperature     = aperature;
    flags.weak_form     = weak_form;
    flags.lagrange_mult = lagrange_mult;
    flags.cold_plasma   = cold_plasma;   
    
    %simulation properties
    settings.omega      = omega;
    settings.N          = N;
    settings.N_harmmax  = N_harmmax;
    settings.x_wall     = x_wall;
    settings.x_ant      = x_ant;
    settings.N_high_res = N_high_res;
    settings.Nk         = Nk;   %number of k-space sample points
    settings.Np         = Np;   %polynomial order
    settings.zeta       = zeta; %k_perp*rho_L number for throughout domain
    settings.ap         = ap;
    settings.r0         = r0;
end

