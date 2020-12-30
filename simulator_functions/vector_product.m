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
%          Dirk van Eester (d.van.eester@fz-juelich.de)                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [matrix] = vector_product(array1,array2)
%VECTOR_PRODUCT calculates N column vectors (N,3,1)* N row vector (N,1,3) =
% N matrices (N,3,3)
%   Input:  vector 1: Nx3(interpreted as N column vectors)
%   Input:  vector 2: Nx3(interpreted as N row    vectors)
%   Output: Matrices: Nx3x3, with each matrix the product of 
    if(size(array1) ~= size(array2))
        error('vectors are not of same size');
    end

    N = size(array1,1);    
    
    matrix = zeros(N,3,3);
    matrix(:,1,1:3) = array1(:,1).*array2;
    matrix(:,2,1:3) = array1(:,2).*array2;
    matrix(:,3,1:3) = array1(:,3).*array2;
    
end

