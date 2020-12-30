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

function [r_high_res,E,E_high_res,timings] = hp_calc_E_all_orders(profiles,flags,settings)
%CP_CALC_ELECTRIC_FIELD 
%   Detailed explanation goes here
   
    %load constants
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
    x_w          = settings.x_wall;             %wall location
    x_a          = settings.x_ant;              %antenna location
    N_high_res   = settings.N_high_res;         %number of high-resolution points for rendering E-field  
    %extract simulation flags
    curve_k_parr = flags.curve_k_parr;          %curved (1) or straight (0) magnetic field lines 
    aperature    = flags.aperature;             %use aperature-type excitation (1) or current sheet (0)
    weak_form    = flags.weak_form;
    lagrange_mult= flags.lagrange_mult;
    cold_plasma  = flags.cold_plasma;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    de_alias                = 0;        %attempt de-aliasing using 2-3rd rule. Use with lagrange multiplier
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x_L    = min(r);
    x_R    = max(r);
    
    if(mod(N,2)==0)
       error('please use an odd number of modes / gridpoints'); 
    end
    
    Nk     = floor(N/2);
    
    
    L      = x_R-x_L;
    a      = L/2;
    
    h      = L/(N-1);
    L      = L+h;
    r0     = (x_R+x_L)/2;    
    x      = r-x_L;
            
   
    
    k0     = real(omega)/c;
    
    
    
    
   
    
    %old method of acquiring indices for x,y and z entries.
%     indices = 0:3*N-1;
%     indices_x = mod(indices,3)==0;
%     indices_y = mod(indices,3)==1;
%     indices_z = mod(indices,3)==2;

    %new method of acquiring indices for x, y and z entries
    %not necessarily faster, but perhaps more readable
    indices_x               = 1:3:3*N; % 1, 4, 7 ... 3*N-2
    indices_y               = 2:3:3*N; % 2, 5, 8 ... 3*N-1
    indices_z               = 3:3:3*N; % 3, 6, 9 ... 3*N
    
    kn            = (-Nk:Nk)';
    kx_test       = (2*pi*kn)/(L);
    fprintf('Loading dielectric tensor( %.2f s )\n',toc);
    
    %load the dielectric tensor
    epsilons = hp_dielectric_tensor_spline(profiles,flags,settings);
    timings.tensor_loaded = toc;
    %swap the 1st and 2nd dimension (effectively transposes each NxN page)
    epsilons = permute(epsilons,[2 1 3 4]);
    fprintf('Done loading dielectric tensor( %.2f s )\n',toc);
    
    Mn = k0^2*epsilons;
    clearvars epsilons;
    
    %compose D
    D  = complex(zeros(N,N,3,3));
    kx = ones(N,1).*kx_test.';
    ky = k(:,2).*ones(1,N);
    kz = k(:,3).*ones(1,N);
   
    D(:,:,1,1) = -(ky.^2 + kz.^2);
    D(:,:,2,2) = -(kx.^2 + kz.^2);
    D(:,:,3,3) = -(kx.^2 + ky.^2);
    
    D(:,:,1,2) = kx.*ky;
    D(:,:,1,3) = kx.*kz;

    D(:,:,2,1) = ky.*kx;
    D(:,:,2,3) = ky.*kz;

    D(:,:,3,1) = kz.*kx;
    D(:,:,3,2) = kz.*ky;


    if(curve_k_parr)
    %assumes curvature
        D(:,:,1,3) = D(:,:,1,3) - 1i*kz./r;
        D(:,:,2,1) = D(:,:,2,1) - 1i*ky./r;
        D(:,:,2,2) = D(:,:,2,2) + 1i*kx./r;
        D(:,:,3,1) = D(:,:,3,1) + 1i*kz./r;
        D(:,:,3,3) = D(:,:,3,3) + 1i*kx./r - 1./r.^2;
    end
    clearvars kx ky kz;
    fprintf('Calculated curl-curl( %.2f s )\n',toc);
    
    
    
    Mn     = Mn + D;
    clearvars D;
    fprintf('Calculated M( %.2f s )\n',toc);
      Mn     = Mn.*exp(1i*(kx_test+kx_test(N))*x.');
     Mn     = h*fft(Mn);
   
    
    
    %create the system matrix
    if(lagrange_mult)
        matrix = complex(zeros(3*N+4,3*N+4));
    else
        matrix = complex(zeros(3*N  ,3*N  ));
    end
    
    %insert the elements in the system matrix
    for dim1 = 1:3
        for dim2 = 1:3
            matrix(dim1:3:3*N,dim2:3:3*N) = Mn(:,:,dim1,dim2);
        end
     end
%     
%     
%     [KK1,KK] = meshgrid(kx_test,kx_test);
%     figure(9)
%     subplot(1,2,1)
%     surf(KK1,KK,abs(squeeze(Mn(:,:,1,1))).','linestyle','none');
%     title('|Mxx|');
%     xlabel('kx');
%     ylabel('kx''');
%     view(2)
%     subplot(1,2,2)
%     surf(KK1,KK,abs(squeeze(Mn(:,:,1,3))).','linestyle','none');
%     title('|Mxz|');
%     xlabel('kx');
%     ylabel('kx''');
%     view(2);
%     DIA
    clearvars Mn;
    fprintf('Created Matrix ( %.2f s )\n',toc);
        
    if(aperature)
        rhs   = zeros(N*3,1);
        E_y_L = 0;
        E_y_R = 1;
        E_z_L = 0;
        E_z_R = 0;
    else
        if(weak_form)
            rhs_r = -exp(-1i*kx_test*(r0+x_a-x_L))*1i*real(omega)*mu0*[0 1 0];
        else
            rhs_r = -1i*real(omega)*mu0*Jant;
        end 
        rhs   = reshape(rhs_r.',[],1);      
        E_y_L = 0;
        E_y_R = 0;
        E_z_L = 0;
        E_z_R = 0;        
    end
    
    
    
    
    
    
    %Set Boundary Conditions
    [~,ind_x_w]    = min(abs(r-(r0+x_w)));
    ind_wall_right = ind_x_w;
    d_ind_wall     = N-ind_x_w+1;
    ind_wall_left  = d_ind_wall;
    
    x_w_L          = x(ind_wall_left );
    x_w_R          = x(ind_wall_right);

    Aeq      = zeros(3*N,4); %the boundary condition matrix
    beq      = zeros(4,1);   %the boundary condition values
    ind_bndc = zeros(4,1);
    
    %Ey at left "wall"   
    Aeq(indices_y,1) = exp(1i*kx_test.*x_w_L);
    beq(1)           = E_y_L;
    ind_bndc(1)      = 3*(ind_wall_left -1)+2;
    
    %Ey at right "wall"
    Aeq(indices_y,2) = exp(1i*kx_test.*x_w_R);
    beq(2)           = E_y_R;
    ind_bndc(2)      = 3*(ind_wall_right-1)+2;
    
    %Ez at left "wall"
    Aeq(indices_z,3) = exp(1i*kx_test.*x_w_L);
    beq(3)           = E_z_L;
    ind_bndc(3)      = 3*(ind_wall_left -1)+3;
    
    %Ez at right "wall"
    Aeq(indices_z,4) = exp(1i*kx_test.*x_w_R);
    beq(4)           = E_z_R;
    ind_bndc(4)      = 3*(ind_wall_right-1)+3;
    
    
    
    
    %apply boundary conditions to the matrix
    if(lagrange_mult)
        %add boundary conditions to the matrix using the Lagrange multiplier
        %method:
    
        %|A   Aeq'||x|   |rhs|
        %|Aeq 0   ||L| = |bnc|
                
        ind_bndc               = 3*N+(1:4);        
        matrix(ind_bndc,1:3*N) = Aeq.';
        matrix(1:3*N,ind_bndc) = Aeq;
        rhs   (ind_bndc)       = beq;        
    else
        matrix(ind_bndc,1:3*N) = Aeq.';
        rhs   (ind_bndc)       = beq;
    end
    
    %do pre-solve anti aliasing (is this the right way to do it?)
    if(de_alias && lagrange_mult)
        matrix_old = matrix;
        K = 1/3*Nk;
        matrix     = zeros(3*N+4+6*K);
        matrix(1:3*N+4,1:3*N+4) = matrix_old;

        nstart = 3*N+4;
        for i = 1:3*K
            matrix(nstart+i,:) = 0;
            matrix(nstart+i,i) = 1;
            matrix(:,nstart+i) = 0;
            matrix(i,nstart+i) = 1;
            rhs   (nstart+i  ) = 0;        
        end

        nstart = 3*N+4+3*K;
        for i = 1:3*K
            matrix(nstart+i,:)       = 0;
            matrix(nstart+i,3*N-i+1) = 1;
            matrix(:,nstart+i)       = 0;
            matrix(3*N-i+1,nstart+i) = 1;
            rhs   (nstart+i  )       = 0;        
        end
    end
    

%     
%      spy(matrix)
%      DIE
%     
    clearvars -except matrix rhs N de_alias N_high_res L kx_test x x_L x_R kn timings;
    fprintf('Matrix system complete  ( %.2f s ), performing matrix inverse\n',toc);
    timings.setup = toc;
    solution = matrix\rhs;    
    timings.solve = toc;
    E_spectrum = solution(1:3*N);   
    
    lambdas    = solution(3*N+1:end);
    fprintf('Matrix inverse complete ( %.2f s ), starting reconstruction of electric field\n',toc);
    
    E_spectrum = reshape(E_spectrum.',[3,N]).';
    
    if(de_alias)
        K = 2/3*Nk;
        E_spectrum(kn<-K,:) = 0;
        E_spectrum(kn> K,:) = 0;
    end
    
    
    
    
    %Definition compatible with MATLAB ifft: xj = (j-1)*L/Nx, j = 1...Nx
    x_high_res = ((1:N_high_res)'-1)*L/(N_high_res);
    
    %old electric field reconstruction method
%     E_high_res = zeros(N_high_res,3);
%     E          = zeros(N,3);   
%     for ind = 1:N
%        kxtn       = kx_test(ind);
%        E          = E         +E_spectrum(ind,:).*exp(1i*x         *kxtn);
%        E_high_res = E_high_res+E_spectrum(ind,:).*exp(1i*x_high_res*kxtn);
%     end
    
    
    %new method based on built-in inverse Fourier transform, much faster
    %for high N_high_res
    E_high_res = N_high_res*ifft(E_spectrum,N_high_res).*exp(1i*kx_test(1)*x_high_res);
    E          = N         *ifft(E_spectrum           ).*exp(1i*kx_test(1)*x         );

    r_high_res      = x_high_res+x_L;
    
    fprintf('Field reconstruction complete ( %.2f s )\n',toc);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%         [rr,kk] = meshgrid(kx_test,x);
%     
%     figure(101)
%     surf(rr,kk,log10(abs(matrix(indices_x,indices_x))),'Linestyle','None')
%     
%     figure(102)
%     surf(rr,kk,abs(matrix(indices_x,indices_y)),'Linestyle','None')
%     figure(103)
%     surf(rr,kk,abs(matrix(indices_x,indices_z)),'Linestyle','None')
%     
%     figure(104)
%     surf(rr,kk,abs(matrix(indices_y,indices_x)),'Linestyle','None')
%     figure(105)
%     surf(rr,kk,abs(matrix(indices_y,indices_y)),'Linestyle','None')
%     figure(106)
%     surf(rr,kk,abs(matrix(indices_y,indices_z)),'Linestyle','None')
%     
%     figure(107)
%     surf(rr,kk,abs(matrix(indices_z,indices_x)),'Linestyle','None')
%     figure(108)
%     surf(rr,kk,abs(matrix(indices_z,indices_y)),'Linestyle','None')
%     figure(109)
%     surf(rr,kk,log10(abs(matrix(indices_z,indices_z))),'Linestyle','None')
%     
%     A = STOP
    
    
    figure(21)
    plot_components_and_abs(kn,log10(abs(E_spectrum)),'-',get(gca,'colororder'));
    xlabel('n');
    ylabel('^{10}log|En|');
    title('Electric field spectrum');
%     figure(22);
%     subplot(2,1,1)
%     plot_components_and_abs(kn,real(E_spectrum),'-',get(gca,'colororder'));
%     subplot(2,1,2)
%     plot_components_and_abs(kn,imag(E_spectrum),'-',get(gca,'colororder'));
    timings.end = toc;
    clearvars -except r_high_res E E_high_res timings;
    
    
end

