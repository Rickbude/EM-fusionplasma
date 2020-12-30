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

clearvars;
close all;

%a very small script that generates input.mat, the file that contains the 
%filename of the input file that should be loaded in input_files
%accepts files with, or without folder specification (default; search in ./
%and in ./input_files).
while(1)
    prompt = 'type the name of the input file you want to load (select a file in ./inputs) \n';
    str = input(prompt,'s');

    input_file_1 = strcat('input_files/',str,'.m');
    input_file_2 = strcat('input_files/',str     );
    input_file_3 = strcat('./',str,'.m');
    input_file_4 = strcat('./',str     );
    %accept both the filename with and without extension
    if(exist(input_file_1,'file'))   
       input_file = input_file_1;
       break;
    elseif(exist(input_file_2,'file'))        
       input_file = input_file_2;
       break;
    elseif(exist(input_file_3,'file'))        
       input_file = input_file_3;
       break;
    elseif(exist(input_file_4,'file'))        
       input_file = input_file_4;
       break;
    else
       fprintf('%s does not exist, please try again \n',input_file); 
    end

end
fprintf('%s is now selected as the default input file\n',input_file);  
save('input_file','input_file');


