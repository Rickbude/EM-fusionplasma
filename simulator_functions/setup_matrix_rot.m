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
function [matrix,rhs_flat] = setup_matrix_rot(A,rhs)
%SETUP_MATRIX Summary of this function goes here
%   Setup the system matrix under the assumption that periodic boundaries
%   can be assumed. A contains the contributions to the matrix

    %extract proper dimensions
    N_coeff = size(A,1);    %number of contributions per equation per row
    N_rows  = size(A,2);    %Number of gridpoints / rows
    N_eq    = size(A,3);    %Number of concurrent equations, e.g. Ex Ey and Ez
    
    %these will hold the indices for the sparse matrix construction
    indexi  = zeros(size(A));   
    indexj  = zeros(size(A));
    
    %column positions relative to main diagonal
    %Assume odd number of coefficients per equation, centered around diagonal
    col_coeff    = -floor(N_coeff/2):floor(N_coeff/2);
    
    %row indices
    ind_row = (1:N_rows) + zeros(N_coeff,N_rows);
    %column indices
    ind_col = (1:N_rows) + ones (N_coeff,N_rows).*col_coeff.';
    
    %correct for out-of-bound contributions
    ind_col(ind_col<1     ) = ind_col(ind_col<1     ) + N_rows ;
    ind_col(ind_col>N_rows) = ind_col(ind_col>N_rows) - N_rows ;
       
    %calculate the actual positions in the system matrix
    %assume the matrix should be "weaved" like:
    %[Ex11 Ey11 Ez12 ...
    % Ex21 Ey21 Ez21 ...
    % .....
    for dim1 = 1:N_eq
        for dim2 = 1:N_eq
            indexi(:,:,dim1,dim2) = dim1 + N_eq*(ind_row-1);
            indexj(:,:,dim1,dim2) = dim2 + N_eq*(ind_col-1);
        end
    end
    
    %reshape into column vectors
    indexi = reshape(indexi,[],1);
    indexj = reshape(indexj,[],1);
    values = reshape(A     ,[],1);
       
    %generate sparse matrix
    matrix = sparse(indexi,indexj,values,N_eq*N_rows,N_eq*N_rows);
    
    %reshape RHS into column vector as well
    rhs_flat = reshape(rhs.',[],1);
    
    clearvars -except rhs_flat matrix;

end

