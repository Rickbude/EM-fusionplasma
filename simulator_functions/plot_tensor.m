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

function [] = plot_tensor_components(r,epsilon)
    
    colors = get(gca,'colororder');
    counter = 1;
    plotnum = 1;
    for dim1 = 1:3
        for dim2 =1:3
            subplot(3,3,plotnum);
            plot(r,epsilon(:,dim1,dim2),'Color',colors(counter,:)), hold on;
            
            counter = counter+1;
            if(counter>length(colors))
                counter = 1;
            end
            plotnum=plotnum+1;
            xlabel('k_x')
            ylabel(sprintf('\\epsilon_{%i%i}',dim1,dim2));
           % ylim([min(epsilon(:,dim1,dim2))-eps max(epsilon(:,dim1,dim2))+eps]);            
        end
    end


end

