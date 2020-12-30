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

function [matrix,rhs_flat] = setup_matrix(columns,rhs,row_indices,column_indices,N)
%SETUP_MATRIX Summary of this function goes here
%   Detailed explanation goes here

%prune away empty columns
empty_rows     = sum(abs(columns),1)==0;
column_indices = column_indices(~empty_rows);
row_indices    = row_indices(~empty_rows);
columns        = columns(:,~empty_rows);


N_col  = length(column_indices);

N_eq   = max(row_indices)+1;

indexi   = zeros(N_col*(N-2),1);
indexj   = zeros(N_col*(N-2),1);

for col_ind = 1:N_col
    %offsets with respect to diagonal
    col_offset = column_indices(col_ind);
    row_offset = row_indices(col_ind);
    
    %the actual row coordinates
    index      = N_eq*(1:N-2)+1+row_offset;
    
    %index required for the "reshape"
    start_ind  = 1+(N-2)*(col_ind-1);
    stop_ind   = start_ind+N-3;
    
    indexi(start_ind:stop_ind) = index;
    indexj(start_ind:stop_ind) = index+col_offset;

end
%values in the matrix
values   = reshape(columns(2:N-1,:),1,[]);
%values in the rhs
rhs_flat = reshape(rhs.',[],1);

%clear out "illegal entries"
valid_values = bitand(bitand(indexi>0,indexj>0),bitand(indexi<3*N+1,indexj<3*N+1));

indexi = indexi(valid_values);
indexj = indexj(valid_values);
values = values(valid_values);

clearvars -except rhs_flat indexi indexj values;
%most time is spent in setting up the sparse matrix



matrix = sparse(indexi,indexj,values);

clearvars -except rhs_flat matrix;

end

