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
function [E] = cp_calc_electric_field_1_5D(r,k,curve_k_parr,omega,B,m,n,q,Jant)
%CP_CALC_ELECTRIC_FIELD 
%   Detailed explanation goes here
    constants;
    
    ky     = k(:,2);
    kz     = k(:,3);    
    
    N      = length(r);
    h      = abs(max(r)-min(r))/(N-1);
    k0     = real(omega)/c;
   
    %epsilon   = cp_dielectric_tensor(omega,B,m,n,q);
    epsilon = cp_dielectric_tensor(omega,B,m,n,q);
    
    %contribitions without derivatives with respect to x
    D0              = zeros(size(epsilon));
    D0(:,1,1)       = -(ky.^2 + kz.^2);
    D0(:,2,2)       = -kz.^2;
    D0(:,3,3)       = -ky.^2;
    D0(:,2,3)       = ky.*kz;
    D0(:,3,2)       = ky.*kz;
    
    %contributions with first derivatives with respect to x
    D1              = zeros(size(epsilon));
    D1(:,1,2)       = -1i*ky;
    D1(:,2,1)       = -1i*ky;
    D1(:,1,3)       = -1i*kz;  
    D1(:,3,1)       = -1i*kz;
    
    %contributions with second derivatives with respect to x
    D2              = zeros(size(epsilon));
    D2(:,2,2)       = 1;
    D2(:,3,3)       = 1;
    
    %extra terms due to curvature
    if(curve_k_parr)
        D0(:,3,1)      = D0(:,1,3) - 1i*kz./r;
        D0(:,1,3)      = D0(:,3,1) + 1i*kz./r;
        D0(:,3,3)      = D0(:,3,3) - 1./(r.^2);
        D0(:,2,1)      = D0(:,2,1) - 1i*ky./r;
        
        D1(:,2,2)      = D1(:,2,2) + 1./r;
        D1(:,3,3)      = D1(:,3,3) + 1./r;
    end
    
    
    %3x3 tensors to be multiplied with E_n-1, E_n en E_n+1
    %A first order finite difference scheme is used
    A_nmin         =  1/h^2*D2 -1/(2*h)*D1;             %E_n-1
    A_n            = -2/h^2*D2 + D0 + k0^2*epsilon;                    %E_n
    A_nplus        =  1/h^2*D2 +1/(2*h)*D1;             %E_n+1

    
    %start building the system matrix
    row_indices     = [0*ones(9,1) 1*ones(9,1) 2*ones(9,1)];    %"equation number"/row (x=0,y=1,z=2)
    column_indices  = [-3:5        -4:4        -5:3       ];    %column relative to Ex_n / Ey_n / Ez_n

    columns = zeros(N,27);
    
    %reshape matrices, suitable for matrix setup
    
    %9 contributions for J_ant_x (eq 1)
    columns(:,1:3  )   =  A_nmin (:,1,1:3);   %E_n-1
    columns(:,4:6  )   =  A_n    (:,1,1:3);   %E_n
    columns(:,7:9  )   =  A_nplus(:,1,1:3);   %E_n+1
    
    %9 contributions for J_ant_y (eq 2)
    columns(:,10:12)   =  A_nmin (:,2,1:3);   %E_n-1
    columns(:,13:15)   =  A_n    (:,2,1:3);   %E_n
    columns(:,16:18)   =  A_nplus(:,2,1:3);   %E_n+1
    
    %9 contributions for J_ant_z (eq 3)
    columns(:,19:21)   =  A_nmin (:,3,1:3);   %E_n-1
    columns(:,22:24)   =  A_n    (:,3,1:3);   %E_n
    columns(:,25:27)   =  A_nplus(:,3,1:3);   %E_n+1

    %rhs: contribution of the antenna current
    rhs          = -1i*real(omega)*mu0*Jant;

    %builds the system matrix and rhs
    [matrix,rhs_column] = setup_matrix(columns,rhs,row_indices,column_indices,N);

    %set the boundary conditions
    
    %Ex (no boundary condition, use backward and forward difference to 
    %"fill the gaps" at the edges. This is not a very rigorous approach
    matrix(1,1:3)               = squeeze(-3/2*D1(1  ,1,1:3)/h + A_n(1,1,1:3) );
    matrix(1,4:6)               = squeeze(   2*D1(2  ,1,1:3)/h);
    matrix(1,7:9)               = squeeze(-1/2*D1(3  ,1,1:3)/h);

    matrix(3*N-2,3*N-8:3*N-6)   = squeeze(1/2*D1(N-2,1,1:3)/h);
    matrix(3*N-2,3*N-5:3*N-3)   = squeeze( -2*D1(N-1,1,1:3)/h);
    matrix(3*N-2,3*N-2:3*N  )   = squeeze(3/2*D1(N  ,1,1:3)/h  + A_n(N,1,1:3) );
    
    rhs_column(1)        = rhs(1,1);
    rhs_column(3*N-2)    = rhs(N,1);   
    
    %other possibility: enforce continuity Ex at the boundary (E1 = E2 and
    %En-1 = En) (third possibility: enforce continuity of dEx?)
    %     matrix(1,1) = 1;
    %     matrix(1,4) = -1;
    %     
    %     matrix(3*N-2,3*N-2)  = 1;
    %     matrix(3*N-2,3*N-5)  = -1;
    
    %Ey = 0 at walls
    matrix(2,2  )       = 1;
    matrix(3*N-1,3*N-1) = 1;

    %Ez = 0 at walls
    matrix(3,3)         = 1;
    matrix(3*N  ,3*N  ) = 1;

    %solve
    E = matrix\rhs_column;
   
    % old method of reshaping
%     indices = 0:3*N-1;
%     Ex = E(mod(indices,3)==0);
%     Ey = E(mod(indices,3)==1);
%     Ez = E(mod(indices,3)==2);
%     E  = [Ex Ey Ez];

    % reshape from [Ex1; Ey1; Ez1; Ex2; Ey2; Ez2;...]    (3Nx1)
    %           to [Ex1  Ey1  Ez1; Ex2  Ey2  Ez; ...]    (Nx3)
    E = reshape(E.',[3,N]).';
end

