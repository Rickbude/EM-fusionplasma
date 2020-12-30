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
function savegoodplot(filename,papersize, margin)
    % function which produces a nice-looking plot
    % and sets up the page for nice printing
    if nargin == 1
        papersize = [15 15];
        margin = 1;
    elseif nargin == 2
        margin = 1;    
    end
    
    height = papersize(1);
    width  = papersize(2);
    
    grid on;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize', [height width]);
    set(gcf,'PaperPosition',[margin margin height-2*margin width-2*margin]);
    set(gcf,'PaperPositionMode','Manual');
    set(gca,'FontSize', 10);
    set(gca,'fontname','times');
    tempfilename = strcat(filename,'.eps');
    print(gcf,'-dpdf', '-r600',tempfilename);
    tempfilename = strcat(filename,'.pdf');
    print(gcf,'-dpdf', '-r600',tempfilename);
end