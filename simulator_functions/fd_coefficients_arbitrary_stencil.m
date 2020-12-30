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

function [coeff] = fd_coefficients_arbitrary_stencil(m,s)
%FD_COEFFICIENTS_ARBITRARY_STENCIL Summary of this function goes here
% Calculate the finite difference coefficients for an arbitrary stencil s
% for an arbitry derivate m
% 
% source: https://en.wikipedia.org/wiki/Finite_difference_coefficient
    ncoeff = length(s);
    %https://en.wikipedia.org/wiki/Finite_difference_coefficient
    %for symmetric stencils (-N/2:N/2), this method is accurate until about
    %m~20
    %convert to vpa for higher derivatives. it will be very slow however 
    %for large stencils!
    
    matrix = zeros(ncoeff);  
%     matrix = vpa(zeros(ncoeff));  
    
    
    for i = 0:ncoeff-1
        for j = 1:ncoeff      
            matrix(i+1,j) = s(j)^(i);
        end
    end        
    
    rhs       = zeros(ncoeff,1);
    rhs(m+1)  = factorial(m);

     
%     rhs       = vpa(zeros(ncoeff,1));
%     rhs(m+1)  = factorial(vpa(m)); 

         
    
    %solve to find the coefficients, convert back to double
    coeff   = double(matrix\rhs);
   
    %upper estimate for the accuracy
    accuracy = ncoeff-m+1;
    
    %"round" the coefficients   
    coeff    = double(int64(coeff*factorial(accuracy)))/factorial(accuracy);
    
            
end

