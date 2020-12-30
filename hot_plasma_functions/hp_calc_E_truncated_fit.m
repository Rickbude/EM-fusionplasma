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

function [E,timings] = hp_calc_E_truncated_fit(profiles,flags,settings)
%CP_CALC_ELECTRIC_FIELD 
%   Detailed explanation goes here

    %% load constants
    constants;  
    
    %% extract r and r-dependent profiles
    B            = profiles.Bb;                 %static magnetic background field profile
    k            = profiles.k;                  %wavenumber (kx, ky, kz)
    m            = profiles.m;                  %mass of all species involved
    n            = profiles.n;                  %density profile of all species involved
    q            = profiles.q;                  %charge of all species involved
    T            = profiles.T;                  %temperature profile of all species involved
    r            = profiles.r;                  %radial position of all gridpoints (relative to the device axis)
    Jant         = profiles.Jant;               %antenna current profile
    
    %% extract simulation settings
    omega        = settings.omega;              %antenna frequency
    N            = settings.N;                  %number of gridpoints
    N_orders     = settings.N_harmmax;          %number of terms in the dielectric tensor sum   
    x_w          = settings.x_wall;             %wall location
    Nk           = settings.Nk;                 %number of k-space sample points
    Np           = settings.Np;                 %polynomial order
    zeta         = settings.zeta;               %k_perp*rho_L number for throughout domain
    
    %% extract simulation flags
    curve_k_parr = flags.curve_k_parr;          %curved (1) or straight (0) magnetic field lines 
    aperature    = flags.aperature;             %use aperature-type excitation (1) or current sheet (0)    
    
    %% debug flags
    %set some flags that can be used to plot various variables in this 
    %function, for debug purposes
    plot_spectrum   = false;
    plot_error      = false;
    plot_coeff      = false;
    plot_M          = false;
    plot_A          = false;
    plot_waitbar    = false;
    
    %Scaling of Ez with respect to the other fields. This scaling is 
    %reversed after the solve phase. facrenorm=1 applies no scaling
    facrenorm       = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    kx     = 1*ones(N,1);
    ky     = k(:,2);
    kz     = k(:,3); 
    
    k0     = real(omega)/c;
    
    timings.start =toc;
    t_elapsed = toc;
    fprintf('Start calculating dielectric tensor (%.0f ms)\n',toc*1000);

   
    
    x_R     = max(r);
    x_L     = min(r);
    r0      = (x_R+x_L)/2;
    L       = x_R-x_L;
    x       = r-r0;
    h       = L/(N-1);
    L       = L+h;
   
    
    %calculate the larmour radius of each species
    rho_L     = larmour_radius(T,m,q,B);
    %find the smallest larmour radius, excluding the electron radius
    rho_L_min = min(rho_L(:,2:end),[],2);
    
    kxrho = zeta./rho_L_min;
    kxlim = pi*N/L*ones(N,1);
    kxabs = 100*ones(N,1);
    %kxmax = kxrho;
    %kxmax = kxabs;
    kxmax = min(kxrho,kxlim);
    
    %kxmax = 500*zeta*ones(N,1);
    
    kx_compressed   = linspace(-1,1,Nk)';
    kx              = kx_compressed;
    
       
    ak = zeros(N,3,3,Np+1); 
    if(plot_waitbar)
        w      = waitbar(0,'loading tensor');
    end
    tprev  = toc;
    
    %load the dielectric tensor components
    epsilons = zeros(N,Nk,3,3);
    for ind = 1:Nk
        kx_test             = kx_compressed(ind)*kxmax;
        k_test              = [kx_test ky kz];
        epsilons(:,ind,:,:) = hp_dielectric_tensor_swanson(N_orders,omega,k_test,B,T,q,n,m); 
        if(toc-tprev > 0.5 && plot_waitbar)
            waitbar(ind/N,w,sprintf('loading tensor (%i / %i)',ind,Nk));
            tprev = toc;
        end
    end
    timings.tensor_loaded = toc;
    if(plot_waitbar)
        close(w)
    end
    fprintf('Start fitting polynomials through tensor components(%.0f ms)\n',toc*1000);
    
    % Set which components are symmetric at kx = 0
    % the symmetric components will have only even polynomial coefficients
    % the antisymmetric components will have only odd polynomial coefficients
    symmetric = [
        1 1 0
        1 1 0
        0 0 1       
    ];
    
    %Split even and odd polynomial coefficients
    p_all  = 0:Np;                     %all polynomial coefficients
    p_even = p_all(mod(0:Np,2)==0);    %the odd  polynomial coefficients
    p_odd  = p_all(mod(0:Np,2)==1);    %the even polynomial coefficients
    
    %construct Vandermonde matrix for both even and odd coefficients
    %Do this only once to save computational time
    kx_matrix_even  = kx_compressed.^p_even;
    kx_matrix_odd   = kx_compressed.^p_odd ;
    kx_matrix       = kx_compressed.^p_all ;
    
    %a decent amount of time is saved by calculating this (pseudo-)inverse
    %only once, and using it later to execute the fits
    kx_matrix_even_inv = pinv(kx_matrix_even).';
    kx_matrix_odd_inv  = pinv(kx_matrix_odd) .';
    
    %"decompression coefficients", used to "decompress" the coefficients
    %after the fit has been executed
    decompression = kxmax.^p_all;          
    
    %now do the fits
    for dim1 = 1:3
        for dim2 = 1:3
            epsilon      = squeeze(epsilons(:,:,dim1,dim2));  
            %holds the compressed coefficients
            p_compressed = zeros(N,Np+1);            
            
            %the actual fit: multiply epsilon by the (pseudo)inverted 
            %Vandermonde matrices to find the least squared error
            
            %Do a fit with a symmetric polynomial if the tensor component is
            %expected to be symmetric
            %do a fit with an odd polynomial if the tensor component is 
            %expected to be antisymmetric
            if(symmetric(dim1,dim2))
                p_compressed(:,p_even+1) = epsilon*kx_matrix_even_inv;
            else
                p_compressed(:,p_odd+1 ) = epsilon*kx_matrix_odd_inv ;
            end             
            %decompress the coefficients
            p = p_compressed./decompression;
            %the polynomial coefficients are stored in ak
            ak(:,dim1,dim2,:) = p;    
        end           
    end
    
    fprintf('Start composing matrix (%.0f ms)\n',toc*1000);  
    
    %plot the coefficients obtained from the fitting procedure
    if(plot_coeff)
        figure(111)         
        
        for dim1 = 1:3
            for dim2 = 1:3
                legendCell = cellstr(num2str(p_all', '%-d'));
                subplot(3,3,(dim1-1)*3+dim2);
                plot(x,log10(abs(squeeze(ak(:,dim1,dim2,:))))),hold on;    
                legend(legendCell);
                xlabel('r [m]');
                ylabel('log_{10}|a_k|');
                title('fit coefficients');
            end
        end        
    end
    
    %plot the fit error
    if(plot_error)
        figure(112)
        for dim1 = 1:3
            for dim2 = 1:3
                actual = squeeze(epsilons(:,:,dim1,dim2));
                fitted = (squeeze(ak(:,dim1,dim2,:)).*decompression)*kx_matrix.';
                errors = sum(abs(actual-fitted),2)/Nk;
                plot(x,errors),hold on;            
            end
        end
        legend('11', '12','13','21','22','23','31','32','33');
        xlabel('r [m]');
        ylabel('error');
        title('dielectric tensor fit error');
    end
    
    %contains the sub-matrices per n-th derivative (M(1) = 0th order
    %derivative, M(2) is 1st order derivative, etc)
    M = zeros(Np+1,N,3,3);    
       
    for m = p_all
        M(m+1,:,:,:) = k0^2*(1i)^(-m)*squeeze(ak(:,:,:,m+1));
    end
    
    [D0,D1,D2] = double_curl_1_5D(N,r,k,curve_k_parr);
    
    %add contributions from the curl-curl tensor
    M(1,:,:,:) = squeeze(M(1,:,:,:)) + D0;
    M(2,:,:,:) = squeeze(M(2,:,:,:)) + D1;
    M(3,:,:,:) = squeeze(M(3,:,:,:)) + D2;
    
    %plot the coefficients obtained from the fitting procedure
    if(plot_M)
        figure(113)            
        for dim1 = 1:3
            for dim2 = 1:3
                legendCell = cellstr(num2str(p_all', '%-d'));
                subplot(3,3,(dim1-1)*3+dim2);
                plot(x,log10(abs(squeeze(M(:,:,dim1,dim2))))),hold on;    
                legend(legendCell);
                xlabel('r [m]');
                ylabel('log_{10}|M_n|');
                title('derivative coefficients');
            end
        end        
    end
    
    clearvars ak D0 D1 D2;
   
    
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
    
    if(plot_A)
        figure(114)        
        ind = -floor(N_coeff/2):floor(N_coeff/2);
        for dim1 = 1:3
            for dim2 = 1:3
                legendCell = cellstr(num2str(ind', '%-d'));
                subplot(3,3,(dim1-1)*3+dim2);
                plot(x,log10(abs(squeeze(A(:,:,dim1,dim2))))),hold on;    
                legend(legendCell);
                xlabel('r [m]');
                ylabel('log_{10}|A_n|');
                title('matrix coefficients');
            end
        end        
    end
     
    %rhs: contribution of the antenna current
    rhs          = -1i*real(omega)*mu0*Jant;

    %builds the system matrix and rhs
    [matrix,rhs_column] = setup_matrix_rot(A,rhs);
        
  
    fprintf('Matrix construction complete, handling boundaries (%.2f s)\n',toc);  
    
    matrix(:,3:3:3*N) = matrix(:,3:3:3*N)./facrenorm;
       
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

