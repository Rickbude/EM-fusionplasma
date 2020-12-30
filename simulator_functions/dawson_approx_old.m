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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Z] = dawson_approx_old(zeta)
%DAWSON_APPROX Approximates the dawson function using series and 
%   maximum error with current settings is 58ppm at x = +/-0.3

abs_zeta = abs(zeta);
part1 = abs_zeta<=0.1;
part2 = bitand(abs_zeta > 0.1, abs_zeta<=10);
part3 = abs_zeta>10;
% part3 = abs_zeta >6    , abs_zeta<=25);
% part4 = bitand(abs_zeta >25   , abs_zeta<=10000);
% part5 = abs_zeta>10000;


%fprintf('part1: %i, part2: %i, part3: %i, part4: %i, part5: %i\n',nnz(part1),nnz(part2),nnz(part3),nnz(part4),nnz(part5));



zeta_inv = 1./zeta;


zeta1     = zeta    (part1);
zeta2     = zeta    (part2);
zeta_inv3 = zeta_inv(part3);
% zeta_inv4 = zeta_inv(part4);
% zeta_inv5 = zeta_inv(part5);

%squaring seems to be optimized in MATLAB, other (integer) powers not so
%much

%series expansion near the origin
Z1  = zeta1 - 2/3*zeta1.*(zeta1.^2) + 4/15*zeta1.*(zeta1.^2).^2 ;
%Z1 = dawson(zeta1);
%this is already faster than MATLAB's Dawson function
%Z2 = 0.5*sqrt(pi)*exp(-zeta2.^2).*erfi(zeta2);
Z2  = dawson_approx_small_z(zeta2,0.4,75);
%Z2 = dawson(zeta2);
%Z2 = dawson(zeta2);
%asymptotic continuation for medium-large zeta
Z3 = 1/2*zeta_inv3 + 1/4*zeta_inv3.*(zeta_inv3.^2) + 3/8*zeta_inv3.*(zeta_inv3.^2).^2;
%asymptotic continuation for large zeta
%Z(part4) = 1./(2*zeta(part4)) + 1./(4*zeta(part4).^3) ;
% Z4 = 1/2*zeta_inv4 + 1/4*zeta_inv4.*(zeta_inv4.^2);
%Z(part4) = 1./(2*zeta(part4)) + 1./(4*zeta(part4).^3) + 3./(8*zeta(part4).^5) ;
%asymptotic continuation for very large zeta
% Z5 = 1/2*zeta_inv5;

Z        = zeros(size(zeta));
Z(part1) = Z1;
Z(part2) = Z2;
Z(part3) = Z3;
% Z(part4) = Z4;
% Z(part5) = Z5;


end

