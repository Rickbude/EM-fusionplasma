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
function [k] = cp_wavenumber_k_parr(omega,k_parr,B,m,n,q)
%DISPERSION_RELATION calculates the wavenumber 
%   calculates the wavenumber at each position in the plasma, as if
%   the plasma was uniform, based on k_parr

    %constants
    c       = 299792458;        %speed of ligt, m/s
    %vacuum wavenumber
    k0      = real(omega)/c;
        
    %load Stix components
    [~,~,P,S,D]=cp_Stix_components(omega,B,m,n,q);
    
    n_parr = k_parr./k0;
    %solve dispersion relation according to Swanson
    
    % As n_perp^4 + Bs n_perp ^2 + Cs = 0
    As = S;
    Bs = D.^2 - S.^2 - P.*S + (P+S).*n_parr.^2;
    Cs = P.*(S.^2 - D.^2 - 2*S.*n_parr.^2 + n_parr.^4);
    
    %normalized wavenumber
    %n21 = (B-sqrt(B.^2-4*A.*C))./(2*A);
    %n22 = (B+sqrt(B.^2-4*A.*C))./(2*A);
    %n2 = real(n21 + n22);
    
    Fs_sq      = Bs.^2-4*As.*Cs;
    
    
   
    n_sq1   = (-Bs-sqrt(Fs_sq))./(2*As);
    n_sq2   = (-Bs+sqrt(Fs_sq))./(2*As);
    
    n1      = sqrt(n_sq1);
    n2      = sqrt(n_sq2);
    n       = sort([n1 n2],2);
    
    
    k       = n.*k0;
        
end

