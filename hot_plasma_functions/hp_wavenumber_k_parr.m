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

function [k_perp] = hp_wavenumber_k_parr(profiles,flags,settings)
%DISPERSION_RELATION calculates the wavenumber 
%   calculates the wavenumber at each position in the plasma, as if
%   the plasma was uniform, based on k_parr

    omega        = settings.omega;
    B            = profiles.Bb;
    k            = profiles.k;
    m            = profiles.m;
    n            = profiles.n;
    q            = profiles.q;
    T            = profiles.T;
    r            = profiles.r;
    %constants
    c       = 299792458;        %speed of ligt, m/s
    %vacuum wavenumber
    k0      = real(omega)/c;
    Nkr     = 500;
    Nki     = 50;
    N       = settings.N;
        
    kx_test_real = [fliplr(-logspace(0,4,Nkr)) logspace(0,4,Nkr)];
    kx_test_imag = [fliplr(-logspace(0,4,Nki)) logspace(0,4,Nki)];
   % kx_test_imag = 0;
    
    [kxr,kxi] = meshgrid(kx_test_real,kx_test_imag);
    kxa = kxr + 1i*kxi;
    
    kxt = reshape(kxa,[],1);
    %kxt = kx_test_real;
    profiles_mod = profiles;
    
    determinants = zeros(N,length(kxt));
    
    w = waitbar(0,'load');
    i_tot = length(kxt);
    for i = 1:N
        waitbar(i/N,w,sprintf('calc dispersion %i/%i',i,N));
        kx = kxt;
        ky = k(i,2)*ones(i_tot,1);
        kz = k(i,3)*ones(i_tot,1);
        kt = [kx ky kz];
        %size(k)
        epsilon = hp_dielectric_tensor_swanson(4,omega,kt,B(i,:),T(i,:),q(i,:),n(i,:),m(i,:));
        
        K        = zeros(i_tot,3,3);
    
        K(:,1,1) = -(ky.^2 + kz.^2);
        K(:,1,2) = kx.*ky;
        K(:,1,3) = kx.*kz;

        K(:,2,1) = ky.*kx;
        K(:,2,2) = -(kx.^2 + kz.^2);
        K(:,2,3) = ky.*kx;

        K(:,3,1) = kz.*kx;
        K(:,3,2) = kz.*kx;
        K(:,3,3) = -(kx.^2 + ky.^2);   
    
    
        M = K + k0^2*epsilon;
        
        %hardcoded 3x3 determinant
        determinants(i,:) = M(:,1,1).*M(:,2,2).*M(:,3,3) ...
                          + M(:,1,2).*M(:,2,3).*M(:,3,1) ...
                          + M(:,1,3).*M(:,2,1).*M(:,3,2) ...
                          - M(:,1,3).*M(:,2,2).*M(:,3,1) ...
                          - M(:,1,2).*M(:,2,1).*M(:,3,3) ...
                          - M(:,1,1).*M(:,2,3).*M(:,3,2);
%         for j = 1:i_tot
%             determinants(i,j) = det(squeeze(matrices(j,:,:)));            
%         end
        
    end
    close(w);
    determinants = reshape(determinants,[N length(kx_test_real) length(kx_test_imag)]);
    
%     figure(3)
%     all_det = squeeze(real(determinants(100,:,:)));
%     pos_det = abs(all_det);
%     neg_det = abs(all_det);
%     pos_det(all_det<0)=0;
%     neg_det(all_det>0)=0;
   
%     size(squeeze(real(determinants(100,:,:))))
%     surf(log10(pos_det))
%     
%     figure(5)
%     size(squeeze(real(determinants(100,:,:))))
%     surf(log10(neg_det))
   
%%%%%%%%%%%%%5
nodes  = zeros(1,4);
solutions = zeros(size(determinants));

k_perp = zeros(N,4*Nkr*Nki);
max_counter = 1;
for n = 1:N
    counter = 1;
    for i = 2:2*Nkr-1         
        for j = 2:2*Nki-1
            node = determinants(n,i,j);
            nodes(1) = determinants(n,i+1,j  );
            nodes(2) = determinants(n,i-1,j  );
            nodes(3) = determinants(n,i,  j-1);
            nodes(4) = determinants(n,i,  j+1);
            %there is a local minimum
            all_nodes = [nodes node];
            if(nnz(abs(nodes)>abs(node)) == 4) 
                %the minimum is at a transition point
                
                
                
                if(~(nnz(real(nodes)==real(node))==4))  
                    
                    solutions(n,i,j)  = 1;
                    k_perp(n,counter) = kxa(j,i);
                    counter = counter + 1;
                    if(counter>max_counter)
                        max_counter = counter;
                    end
                end
            end                    
        end
    end
end
disp(max_counter)
k_perp = k_perp(:,1:max_counter);

% disp(size(determinants));
% size(solutions);
% figure(1)
% surf(log10(abs(squeeze(determinants(400,:,:)))));
% figure(2)
% spy(squeeze(solutions(400,:,:)));
% 
% DIE;

%%%%%%%%%%%%%%%%%%%%%%%
% 
%     k_perp = zeros(N,i_tot);
%     for i = 1:N
%         crossings1 = bitand(sign(real(determinants(i,2:end))) == -sign(real(determinants(i,1:end-1))),sign(real(determinants(i,2:end)))~=0);
%        % crossings2 = bitand(sign(imag(determinants(i,2:end))) == -sign(imag(determinants(i,1:end-1))),sign(imag(determinants(i,2:end)))~=0);
% 
%        % ks(counter,1:(nnz(crossings1)+nnz(crossings2))) = [kx(crossings1); kx(crossings2)]';
% 
%         k_perp(i,1:(nnz(crossings1))) = kxt(crossings1);
%        % k_perp(counter,1:(nnz(crossings2))) = kx(crossings2);
%     
%     end

end

