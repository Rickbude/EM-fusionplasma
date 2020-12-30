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

function [] = plot_all_tensor_components(r,epsilon,d_epsilon,H_epsilon)
%PLOT_ALL_TENSOR_COMPONENTS Summary of this function goes here
%   Detailed explanation goes here
    f = figure(11);
    f.Name='epsilon'; 
    plot_tensor(r,real(epsilon)),hold on;
    
    f = figure(21);
    f.Name='d/dx epsilon'; 
    plot_tensor(r,real(squeeze(d_epsilon(:,:,:,1)))),hold on;
    f = figure(22);
    f.Name='d/dy epsilon'; 
    plot_tensor(r,real(squeeze(d_epsilon(:,:,:,2)))),hold on;
    f = figure(23);
    f.Name='d/dz epsilon'; 
    plot_tensor(r,real(squeeze(d_epsilon(:,:,:,3)))),hold on;
    
    f = figure(31);
    f.Name='d/dxx epsilon';
    plot_tensor(r,real(squeeze(H_epsilon(:,:,:,1,1))));
    f = figure(32);
    f.Name='d/dxy epsilon';
    plot_tensor(r,real(squeeze(H_epsilon(:,:,:,1,2))));
    f = figure(33);
    f.Name='d/dxz epsilon';
    plot_tensor(r,real(squeeze(H_epsilon(:,:,:,1,3))));
    f = figure(34);
    f.Name='d/dyx epsilon';
    plot_tensor(r,real(squeeze(H_epsilon(:,:,:,2,1))));
    f = figure(35);
    f.Name='d/dyy epsilon';
    plot_tensor(r,real(squeeze(H_epsilon(:,:,:,2,2))));
    f = figure(36);
    f.Name='d/dyz epsilon';
    plot_tensor(r,real(squeeze(H_epsilon(:,:,:,2,3))));
    f = figure(37);
    f.Name='d/dzx epsilon';
    plot_tensor(r,real(squeeze(H_epsilon(:,:,:,3,1))));
    f = figure(38);
    f.Name='d/dzy epsilon';
    plot_tensor(r,real(squeeze(H_epsilon(:,:,:,3,2))));
    f = figure(39);
    f.Name='d/dzz epsilon';
    plot_tensor(r,real(squeeze(H_epsilon(:,:,:,3,3))));
end

