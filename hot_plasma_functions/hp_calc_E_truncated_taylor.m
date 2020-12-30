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

function [E,timings] = hp_calc_E_truncated_taylor(profiles,flags,settings)
%CP_CALC_ELECTRIC_FIELD 
%   Detailed explanation goes here
    constants;
    timings.start =toc;
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
    
    plot_spectrum   = false;
    plot_error      = false;
    plot_coeff      = false;
    plot_M          = false;
    plot_A          = false;
    facrenorm       = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dk     = zeta;
    kx     = 0*ones(N,1);
    ky     = k(:,2);
    kz     = k(:,3); 
    
    k_test = [kx ky kz];    
    
    h      = abs(max(r)-min(r))/(N-1);
    
    k0     = real(omega)/c;
    
    
    t_elapsed = toc;
    fprintf('Start calculating dielectric tensor (%.0f ms)\n',toc*1000);

    
    if(mod(N,2)==0)
       error('please use an odd number of modes / gridpoints'); 
    end
    
    
    x_R     = max(r);
    x_L     = min(r);
    r0      = (x_R+x_L)/2;
    L       = x_R-x_L;
    x       = r-r0;
    h       = L/(N-1);
    L       = L+h;
    [epsilon, d_epsilon, H_epsilon] = hp_dielectric_tensor_and_derivatives(N_orders,omega,k_test,dk,B,T,q,n,m);
    timings.tensor_loaded = toc;
    %plot_all_tensor_components(r,epsilon,d_epsilon,H_epsilon);
    
    
    t_elapsed = toc;
    fprintf('Done calculating dielectric tensor and derivatives (%.0f ms)\n',t_elapsed*1000);

    %load the components that are related to derivatives with respect to x
    [D0,D1,D2] = double_curl_1_5D(N,r,k,curve_k_parr);
        
    scriptE0 = epsilon - kx.*squeeze(d_epsilon(:,:,:,1)) + 1/2.*kx.^2.*squeeze(H_epsilon(:,:,:,1,1)); 
    scriptE1 = squeeze(d_epsilon(:,:,:,1)) - kx.*squeeze(H_epsilon(:,:,:,1,1));
    scriptE2 = 0.5*squeeze(H_epsilon(:,:,:,1,1));
    
    Np = 2;
    
    M  = zeros(Np+1,N,3,3);
    M(1,:,:,:) = D0 + 1 * k0^2*scriptE0;
    M(2,:,:,:) = D1 - 1i* k0^2*scriptE1;
    M(3,:,:,:) = D2 - 1 * k0^2*scriptE2;
        
    %load the central difference coefficients
    [~,coeff_cent,~] = fd_coefficients(Np);
    
    %calculate the contributions for each gridpoint relative to the center
    %gridpoint. 
    % A(1) corresponds to E_-N/2
    % A(N) corresponds to E_+N/2
    N_coeff = size(coeff_cent,2);
    A = zeros(N_coeff,N,3,3);
    for m = 0:Np  %m-th order derivative
        for i = 1:N_coeff %i-th gridpoint (E_n-N ... E_n-1 E_n E_n+1 ... E_N)
            A(i,:,:,:) = A(i,:,:,:) + coeff_cent(m+1,i)*M(m+1,:,:,:)/h^m;
        end
    end
    %A = A./multi.';
    
    
     
    %rhs: contribution of the antenna current
    rhs          = -1i*real(omega)*mu0*Jant;

    %builds the system matrix and rhs
    [matrix,rhs_column] = setup_matrix_rot(A,rhs);
        
  
    fprintf('Matrix construction complete, handling boundaries (%.2f s)\n',toc);  
    
%     figure(99)
%     spy(matrix)
    
    
    clearvars columns row_indices column_indices;

    %%%%%%%%%%%$$$$%%%%%%%% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %Set Boundary Conditions
    [~,ind_x_w]    = min(abs(r-(r0+x_w)));
    ind_wall_right = ind_x_w;
    d_ind_wall     = N-ind_x_w+1;
    ind_wall_left  = d_ind_wall;
    
        
    if(aperature)
        rhs_column = zeros(N*3+1,1);
        E_y_L = 0;
        E_y_R = 1;
        E_z_L = 0;
        E_z_R = 0;
    else       
        E_y_L = 0;
        E_y_R = 0;
        E_z_L = 0;
        E_z_R = 0;        
    end
    
       
    %set the boundary conditions using a Lagrange multiplier type of set-up  
    N_bc_eq = 7;
    Aeq = sparse(N_bc_eq,3*N);
    %Metallic wall / aperature electric field restrictions
    Aeq(1,3*(ind_wall_left -1)+2) = 1;  %Ey left wall
    Aeq(2,3*(ind_wall_left -1)+3) = 1;  %Ez left wall
    Aeq(3,3*(ind_wall_right-1)+2) = 1;  %Ey right wall
    Aeq(4,3*(ind_wall_right-1)+3) = 1;  %Ez right wall
    %periodic boundary conditions
    Aeq(5 ,1  )                   = 1;  %Ex left boundary
    Aeq(6 ,2  )                   = 1;  %Ey left boundary
    Aeq(7 ,3  )                   = 1;  %Ez left boundary
    Aeq(5 ,3*N-2)                 = -1; %Ex right boundary
    Aeq(6 ,3*N-1)                 = -1; %Ey right boundary
    Aeq(7 ,3*N-0)                 = -1; %Ez right boundary
    
    %add boundary conditions to system matrix
    matrix(3*N+(1:N_bc_eq),1:3*N) = Aeq;
    matrix(1:3*N,3*N+(1:N_bc_eq)) = Aeq.';
    
    %add boundary conditions to RHS      
    rhs_column(3*N+1)   = E_y_L;    %Ey left wall
    rhs_column(3*N+2)   = E_z_L;    %Ez left wall
    rhs_column(3*N+3)   = E_y_R;    %Ey right wall
    rhs_column(3*N+4)   = E_z_R;    %Ez right wall
    rhs_column(3*N+5)   = 0;        %Ex periodic
    rhs_column(3*N+6)   = 0;        %Ey periodic
    rhs_column(3*N+7)   = 0;        %Ez periodic
    
    fprintf('Set up matrix, starting solve (%.2f s)\n',toc);
    
    fprintf('fillfactor = %.2f%%\n',nnz(matrix)/N^2);
    
    
    clearvars -except r kx matrix rhs_column N facrenorm plot_spectrum multi timings;
    
    fprintf('Reciprocal matrix condition number %.2e\n',1/condest(matrix));
    
    timings.setup = toc;
    %solve for the electric field    
    E = matrix\rhs_column;
    timings.solve = toc;
    fprintf('Matrix inversion complete (%.2f s)\n',toc); 
    
    %1:3N contains E, 3N+1:3N+N_eq contains the Lagrange multipliers
    E = E(1:3*N);

    % reshape from [Ex1; Ey1; Ez1; Ex2; Ey2; Ez2;...]    (3Nx1)
    %           to [Ex1  Ey1  Ez1; Ex2  Ey2  Ez; ...]    (Nx3)
    E  = reshape(E.',[3,N]).';
    
    
    E(:,3) = E(:,3)./facrenorm;
    
    if(plot_spectrum)
        %plot the spectrum
        Fs         = N/abs(max(r)-min(r))*(2*pi);    %samples per meter
        k_spectrum = linspace(-Fs/2,Fs/2,N)';

        E_spectrum = fftshift(fft(E),1)/N;
        figure(123)
        plot(k_spectrum,abs(E_spectrum));
        xlabel('k');
        ylabel('E_k');
        title('Electric field spectrum');
    end
    timings.end = toc;
    clearvars -except E timings;
end

