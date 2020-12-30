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

%load a "timings" file. Generate it first using compare_3_models.m
load('results/times/500points_fixed_report.mat');
matrix_time = elapsed_solve-elapsed_setup;
other_time  = elapsed_end - matrix_time;
all_time    = elapsed_end - elapsed_start;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%          plot RRSE as function of N      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
loglog(Ns,total_errors)
xlabel('$N$','interpreter','latex')
ylabel('RRSE','interpreter','latex')
title('RRSE as function of $N$','interpreter','latex');
legend('Trunc. fit','Trunc. Taylor','All-orders');
ylim([1e-4 1]);
xlim([100 10000]);
savegoodplot('./plot_out/RRSE_afo_N',[10 8],0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% plot computational time as function of N %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)

colors = get(gca,'colororder');

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
xlabel('$N$','interpreter','latex')
ylabel('$t$ [s]','interpreter','latex')
ylim([0.001 1000]);
grid on;
title('Matrix setup','interpreter','latex');
%legend(sprintf('All-orders $(~N^{%.1f})$',p3(1)),sprintf('polynomial $(~N^{%.1f})$',p1(1)),sprintf('Taylor $(~N^{%.1f})$',p2(1)),'Location','southeast','interpreter','latex');
legend('All-orders','Polynomial','Taylor', 'Location','southwest');

%Add scaling with N to graph
text(3e3,5e2 ,sprintf('$N^{%.1f}$',p3(1)),'interpreter','latex');   %all-orders
text(1e4,2e1 ,sprintf('$N^{%.1f}$',p1(1)),'interpreter','latex');   %polynomial
text(2e4,6e-1,sprintf('$N^{%.1f}$',p2(1)),'interpreter','latex');   %Taylor


set(gca,'FontSize', 10);
set(gca,'fontname','times');

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

%Add scaling with N to graph
text(1e4,5e2 ,sprintf('$N^{%.1f}$',p3(1)),'interpreter','latex');   %all-orders
text(6e3,2e0 ,sprintf('$N^{%.1f}$',p1(1)),'interpreter','latex');   %polynomial
text(2e4,6e-2,sprintf('$N^{%.1f}$',p2(1)),'interpreter','latex');   %Taylor

set(gca,'XScale','log')
set(gca,'YScale','log')
hold off;
xlabel('$N$','interpreter','latex')
ylabel('$t$ [s]','interpreter','latex')
ylim([0.001 1000]);
title('Matrix solve','interpreter','latex');
%legend(sprintf('All-orders $(~N^{%.1f})$',p3(1)),sprintf('polynomial $(~N^{%.1f})$',p1(1)),sprintf('Taylor $(~N^{%.1f})$',p2(1)),'Position',[0.4 0.4 0.2 0.2],'interpreter','latex');
grid on;
set(gca,'FontSize', 10);
set(gca,'fontname','times');
% disp(10.^(polyval(p3,log10([1e3 1e6 1e9]))))
% disp(10.^(polyval(p1,log10([1e3 1e6 1e9]))))
% disp(10.^(polyval(p2,log10([1e3 1e6 1e9]))))


savegoodplot('./plot_out/computational_time_scaling',[10 8],0);





N_est = [2e2 4e4 8e6];

times_est21 = 10.^(polyval(p1,log10(N_est)));
times_est22 = 10.^(polyval(p2,log10(N_est)));
times_est23 = 10.^(polyval(p3,log10(N_est)));

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
