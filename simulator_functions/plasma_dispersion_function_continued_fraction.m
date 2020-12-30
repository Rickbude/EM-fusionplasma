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

function [Z] = plasma_dispersion_function_continued_fraction(zeta,N)
%DAWSON_CONTINUED_FRACTION Calculate the plasma dispersion function for
%large imaginary arguments, using the continued fraction from 
%Fried & Conte 1961: the plasma dispersion function, with the correction
%from 10.2307/2003748 (Review by Y.L.L of the book above)
%   zeta is the argument, N the number of terms
%This works especially well for Im(zeta)>1, for somewhat larger Re(zeta)

    An    = 0;  %A_0   = 0
    Bn    = 1;  %B_0   = 1
    Aprev = 1;  %A_n-1 = 1
    Bprev = 0;  %B_n-1 = 0
    
    zeta2 = zeta.^2;
    
    for n = 0:N
        if(n==0)    
            %a_n+1 = a1 = zeta
            a = zeta;
        else
            %a_n+1 = n(2n-1)/n, n = 1,2,...
            a = n*(2*n-1)/2;
        end
        %b_n+1 = -zeta^2 + 0.5 + 2n, n = 0,1,2,...
        b = -zeta2 + 0.5 + 2*n;
        
        %In the work of Conte:
        %A_n+1 = b_n+1 A_n + a_n+1 A_n-1
        %B_n+1 = b_n+1 A_n + a_n+1 B_n-1
        %Seems like it should be adjusted to:
        Anew = b.*An - a.*Aprev;
        Bnew = b.*Bn - a.*Bprev;
        
        Aprev = An;
        Bprev = Bn;
        
        An = Anew;
        Bn = Bnew;
    end

    %This - is not there in the work of Fried and Conte, but it seems to be
    %necessary
    Z = -An./Bn;
   
    %Calculate results for zeta with negative imaginary part too
    zeta_neg = zeta(imag(zeta)<0);
    if(~isempty(zeta_neg))
        Z(imag(zeta)<0) = 2*1i*sqrt(pi)*exp(-zeta_neg.^2)+plasma_dispersion_function_continued_fraction(conj(zeta_neg),N);
    end
end

