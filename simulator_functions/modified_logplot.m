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
function [outputArg1,outputArg2] = modified_logplot(x,y)
%MODIFIED_LOGPLOT Summary of this function goes here
%   Detailed explanation goes here
    
    ylog = sign(y).*log10(abs(y));    
    yreg = bitand(y>-1,y<1);
    ylog(yreg) = y(yreg);
    
    plot(x,ylog)
    
    if(max(ylog)-min(ylog)<10)
        yticks([-10 -8 -6 -4 -2 0 2 4 6 8 10])
        yticklabels({'-10^{10}','-10^8','-10^6','-10^4','-10^2','0','10^2','10^4','10^6','10^8','10^{10}'})
    
    else
        yticks([-9  -6  -3 0 3 6 9])
        yticklabels({'-10^9','-10^6','-10^3','0','10^3','10^3','10^9'})
    
    end
    


end

