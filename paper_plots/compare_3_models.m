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

addpath('../plasma_property_functions');    % load functions that contain  plasma properties (cyclotron frequency etc)
addpath('../hot_plasma_functions');         % load functions specifically for hot plasma (cp Stix tensor, etc)
addpath('../hot_plasma_functions_clean');   % load functions specifically for hot plasma (cp Stix tensor, etc)
addpath('../cold_plasma_functions');        % load functions specifically for cold plasma (cp Stix tensor, etc)
addpath('../EM_property_functions');        % load functions for EM properties (Poynting vector etc)
addpath('../profile_functions');            % load functions for generating profiles
addpath('../simulator_functions');          % load functions that support running simulations (general)
addpath('../plasma_profiles');              % load parameter files
addpath('../error_functions');              % load different error metrics

clearvars;
clear all;
tic;

load('../input_file');
input_file = strcat('../',input_file);

close all;
profile on;
%load physical constants
constants;

load('results/previous_results_ao_N12601_Nk12601_JET_75H');

N_Ns      = 500;
N_min     = 200;
N_max     = 100000;

%Ns        = double(unique(int64(linspace(N_min,N_max,N_Ns))));
Ns        = double(unique(floor(int64(logspace(log10(N_min),log10(N_max),N_Ns).'/18)*18+1)));
N_Ns      = length(Ns);

total_errors    = zeros(N_Ns,3);
elapsed_start           = zeros(N_Ns,3);
elapsed_tensor_loaded   = zeros(N_Ns,3);
elapsed_setup           = zeros(N_Ns,3);
elapsed_solve           = zeros(N_Ns,3);
elapsed_end             = zeros(N_Ns,3);
w = waitbar(0,'starting calculations');

%fit
counter = 1;
for ind_N = 1:length(Ns)        
        N    = Ns(ind_N);
        waitbar(counter/N_Ns,w,sprintf('Truncated fit: N = %i/%i',N,max(Ns)));        
        override.N    = N;
        override.Nk   = 50;
        [profiles,flags,settings] = load_profile_override(input_file,override);
        r  = profiles.r;
        r0 = (max(r)+min(r))/2;
        ap = settings.ap;

        E_interp  = interp1(r_prev,E_prev,r,'spline');
        mask_zero = abs((r-r0))>=ap;

        settings.N = Ns(ind_N);
        tic;
        [E,timings] = hp_calc_E_truncated_fit(profiles,flags,settings);
        elapsed_solve(ind_N,1) = timings.solve - timings.setup;
        
        model = 1;
        elapsed_start        (ind_N,model) = timings.start;
        elapsed_tensor_loaded(ind_N,model) = timings.tensor_loaded;
        elapsed_setup        (ind_N,model) = timings.setup;
        elapsed_solve        (ind_N,model) = timings.solve;
        elapsed_end          (ind_N,model) = timings.end;
        
        total_error = RRSE(E_interp(~mask_zero,:),E(~mask_zero,:));
        total_errors(ind_N,1) = mean(total_error);
        counter = counter+1;   
end
close(w)
w = waitbar(0,'starting calculations');
%Taylor
counter = 1;
for ind_N = 1:length(Ns)        
        N    = Ns(ind_N);
        waitbar(counter/N_Ns,w,sprintf('Truncated Taylor: N = %i/%i',N,max(Ns)));           
        override.N    = N;
        [profiles,flags,settings] = load_profile_override(input_file,override);
        r = profiles.r;

        E_interp  = interp1(r_prev,E_prev,r,'spline');
        mask_zero = abs((r-r0))>=ap;

        settings.N = Ns(ind_N);
        tic;
        [E,timings] = hp_calc_E_truncated_taylor(profiles,flags,settings);
        model = 2;
        elapsed_start        (ind_N,model) = timings.start;
        elapsed_tensor_loaded(ind_N,model) = timings.tensor_loaded;
        elapsed_setup        (ind_N,model) = timings.setup;
        elapsed_solve        (ind_N,model) = timings.solve;
        elapsed_end          (ind_N,model) = timings.end;
        total_error = RRSE(E_interp(~mask_zero,:),E(~mask_zero,:));
        total_errors(ind_N,2) = mean(total_error);
        counter = counter+1;   
end
close(w)
w = waitbar(0,'starting calculations');
%all-orders
counter = 1;
for ind_N = 1:length(Ns)      
        N    = Ns(ind_N);
        
        %be nice to your PC. 
        % for 16BG RAM, limit to N=5000
        % for 32GB RAM, limit to N=9000
        if(N>9000)
           continue; 
        end
         waitbar(counter/N_Ns,w,sprintf('all-orders: N = %i/%i',N,max(Ns)));        
        override.N    = N;
        override.Nk   = max(N_min,ceil(N/10));    %
        [profiles,flags,settings] = load_profile_override(input_file,override);
        r = profiles.r;

        E_interp  = interp1(r_prev,E_prev,r,'spline');
        mask_zero = abs((r-r0))>=ap;

        settings.N = Ns(ind_N);
        tic;
        [r_high_res,E,E_high_res,timings] = hp_calc_E_all_orders      (profiles,flags,settings);
        model = 3;
        elapsed_start        (ind_N,model) = timings.start;
        elapsed_tensor_loaded(ind_N,model) = timings.tensor_loaded;
        elapsed_setup        (ind_N,model) = timings.setup;
        elapsed_solve        (ind_N,model) = timings.solve;
        elapsed_end          (ind_N,model) = timings.end;
        total_error = RRSE(E_interp(~mask_zero,:),E(~mask_zero,:));
        total_errors(ind_N,3) = mean(total_error);
        counter = counter+1;   
end
close(w)

figure(1)
loglog(Ns,total_errors)
xlabel('N')
ylabel('RRSE')
title('RRSE as function of N');
legend('Trunc. fit','Trunc. Taylor','All-orders');
ylim([1e-4 1]);
xlim([100 10000]);
savegoodplot('./plot_out/RRSE_afo_N',[10 8],0);


matrix_time = elapsed_solve-elapsed_setup;
other_time  = elapsed_end - matrix_time;
all_time    = elapsed_end - elapsed_start;

figure(2)

colors = get(gca,'colororder');
N_est = [2e2 4e4 8e6];
subplot(1,2,2)
hold on;
curr_time = matrix_time;
to_fit    = bitand(curr_time~=0,Ns>3000);
p1 = polyfit(log10(Ns(to_fit(:,1))),log10(curr_time(to_fit(:,1),1)),1);
p2 = polyfit(log10(Ns(to_fit(:,2))),log10(curr_time(to_fit(:,2),2)),1);
p3 = polyfit(log10(Ns(to_fit(:,3))),log10(curr_time(to_fit(:,3),3)),1);
plot(Ns,10.^(polyval(p3,log10(Ns))),'Color',colors(3,:));
plot(Ns,10.^(polyval(p1,log10(Ns))),'Color',colors(1,:));
plot(Ns,10.^(polyval(p2,log10(Ns))),'Color',colors(2,:));
plot(Ns,curr_time(:,3),'.','Color',colors(3,:))
plot(Ns,curr_time(:,1),'.','Color',colors(1,:))
plot(Ns,curr_time(:,2),'.','Color',colors(2,:))
set(gca,'XScale','log')
set(gca,'YScale','log')
hold off;
xlabel('N')
ylabel('t [s]')
ylim([0.001 1000]);
title('Time spent in matrix solve');
legend(sprintf('All-orders (~N^{%.1f})',p3(1)),sprintf('polynomial (~N^{%.1f})',p1(1)),sprintf('Taylor (~N^{%.1f})',p2(1)),'Location','southeast');
grid on;
disp('time spent on matrix invert: ');
% disp(10.^(polyval(p3,log10([1e3 1e6 1e9]))))
% disp(10.^(polyval(p1,log10([1e3 1e6 1e9]))))
% disp(10.^(polyval(p2,log10([1e3 1e6 1e9]))))


times_est21 = 10.^(polyval(p1,log10(N_est)));
times_est22 = 10.^(polyval(p2,log10(N_est)));
times_est23 = 10.^(polyval(p3,log10(N_est)));


subplot(1,2,1)
hold on;
curr_time = other_time;
to_fit    = bitand(curr_time~=0,Ns>3000);
p1 = polyfit(log10(Ns(to_fit(:,1))),log10(curr_time(to_fit(:,1),1)),1);
p2 = polyfit(log10(Ns(to_fit(:,2))),log10(curr_time(to_fit(:,2),2)),1);
p3 = polyfit(log10(Ns(to_fit(:,3))),log10(curr_time(to_fit(:,3),3)),1);
plot(Ns,10.^(polyval(p3,log10(Ns))),'Color',colors(3,:));
plot(Ns,10.^(polyval(p1,log10(Ns))),'Color',colors(1,:));
plot(Ns,10.^(polyval(p2,log10(Ns))),'Color',colors(2,:));
plot(Ns,curr_time(:,3),'.','Color',colors(3,:))
plot(Ns,curr_time(:,1),'.','Color',colors(1,:))
plot(Ns,curr_time(:,2),'.','Color',colors(2,:))
set(gca,'XScale','log')
set(gca,'YScale','log')
hold off;
xlabel('N')
ylabel('t [s]')
ylim([0.001 1000]);
grid on;
title('Time spent setting up matrix system');
legend(sprintf('All-orders (~N^{%.1f})',p3(1)),sprintf('polynomial (~N^{%.1f})',p1(1)),sprintf('Taylor (~N^{%.1f})',p2(1)),'Location','southeast');
%legend('All-orders','Polynomial','Taylor');
disp('time spent setting up: ');

times_est11 = 10.^(polyval(p1,log10(N_est)));
times_est12 = 10.^(polyval(p2,log10(N_est)));
times_est13 = 10.^(polyval(p3,log10(N_est)));

times_est1  = (times_est11 + times_est21).*[1 9 81];
times_est2  = times_est12 + times_est22.*[1 3 9];
times_est3  = times_est13 + times_est23;
% disp(datestr(10.^(polyval(p3,log10([1e3 1e6 1e9]))),'HH:MM:SS'))
% disp(10.^(polyval(p1,log10([1e3 1e6 1e9]))))
% disp(10.^(polyval(p2,log10([1e3 1e6 1e9]))))

t = datetime('now') + seconds(times_est1);
D = caldiff([datetime('now') t],{'months','weeks','days','Time'});
disp(D);
t = datetime('now') + seconds(times_est2);
D = caldiff([datetime('now') t],{'months','weeks','days','Time'});
disp(D);
t = datetime('now') + seconds(times_est3);
D = caldiff([datetime('now') t],{'years','months','weeks','days','Time'});
disp(D);


% subplot(1,3,2)
% hold on;
% curr_time = elapsed_tensor_loaded;
% p1 = polyfit(log10(Ns),log10(curr_time(:,1)),1);
% p2 = polyfit(log10(Ns),log10(curr_time(:,2)),1);
% p3 = polyfit(log10(Ns),log10(curr_time(:,3)),1);
% 
% plot(Ns,10.^(polyval(p3,log10(Ns))),'Color',colors(3,:));
% plot(Ns,10.^(polyval(p1,log10(Ns))),'Color',colors(1,:));
% plot(Ns,10.^(polyval(p2,log10(Ns))),'Color',colors(2,:));
% plot(Ns,curr_time(:,3),'.','Color',colors(3,:))
% plot(Ns,curr_time(:,1),'.','Color',colors(1,:))
% plot(Ns,curr_time(:,2),'.','Color',colors(2,:))
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% hold off;
% xlabel('N')
% ylabel('t [s]')
% ylim([0.001 100]);
% title('Time spent loading dielectric tensor');
% legend(sprintf('All-orders (~N^{%.1f})',p3(1)),sprintf('polynomial (~N^{%.1f})',p1(1)),sprintf('Taylor (~N^{%.1f})',p2(1)));
% 

savegoodplot('./plot_out/computational_time_scaling',[12 10],0);

figure(3)
hold on;
colors = get(gca,'colororder');
p1 = polyfit(log10(Ns),log10(other_time(:,1)),1);
p2 = polyfit(log10(Ns),log10(other_time(:,2)),1);
p3 = polyfit(log10(Ns),log10(other_time(:,3)),1);
plot(Ns,10.^(polyval(p1,log10(Ns))),'Color',colors(1,:));
plot(Ns,10.^(polyval(p2,log10(Ns))),'Color',colors(2,:));
plot(Ns,10.^(polyval(p3,log10(Ns))),'Color',colors(3,:));
plot(Ns,other_time(:,1),'.','Color',colors(1,:))
plot(Ns,other_time(:,2),'.','Color',colors(2,:))
plot(Ns,other_time(:,3),'.','Color',colors(3,:))
set(gca,'XScale','log')
set(gca,'YScale','log')
hold off;
xlabel('N')
ylabel('t [s]')
title('Time spent in other code');
legend(sprintf('polynomial (~N^{%.1f})',p1(1)),sprintf('Taylor (~N^{%.1f})',p2(1)),sprintf('All-orders (~N^{%.1f})',p3(1)));

figure(4)
hold on;
colors = get(gca,'colororder');
p1 = polyfit(log10(Ns),log10(all_time(:,1)),1);
p2 = polyfit(log10(Ns),log10(all_time(:,2)),1);
p3 = polyfit(log10(Ns),log10(all_time(:,3)),1);
plot(Ns,10.^(polyval(p1,log10(Ns))),'Color',colors(1,:));
plot(Ns,10.^(polyval(p2,log10(Ns))),'Color',colors(2,:));
plot(Ns,10.^(polyval(p3,log10(Ns))),'Color',colors(3,:));
plot(Ns,all_time(:,1),'.','Color',colors(1,:))
plot(Ns,all_time(:,2),'.','Color',colors(2,:))
plot(Ns,all_time(:,3),'.','Color',colors(3,:))
set(gca,'XScale','log')
set(gca,'YScale','log')
hold off;
xlabel('N')
ylabel('t [s]')
title('Time spent in all code combined');
legend(sprintf('polynomial (~N^{%.1f})',p1(1)),sprintf('Taylor (~N^{%.1f})',p2(1)),sprintf('All-orders (~N^{%.1f})',p3(1)));


figure(5)
loglog(Ns,elapsed_tensor_loaded)
xlabel('N')
ylabel('t [s]')
title('Time spent loading tensor');
legend('Trunc. fit','Trunc. Taylor','All-orders');


figure(6)
loglog(Ns,all_time)
xlabel('N')
ylabel('t [s]')
title('Total time spent'),hold on;
%plot(Ns,0.001*Ns.*log(Ns));
hold off
legend('Trunc. fit','Trunc. Taylor','All-orders');


figure(7)
loglog(Ns,all_time)
xlabel('Resolution / Number of gridpoints')
ylabel('Time spent [s]')
title('Comparison of slow, fast and hybrid approach'),hold on;
%plot(Ns,0.001*Ns.*log(Ns));
hold off
legend('Hybrid','Fast','Slow');

%save('results/times/500points','Ns','total_errors','elapsed_start','elapsed_tensor_loaded','elapsed_setup','elapsed_solve','elapsed_end');