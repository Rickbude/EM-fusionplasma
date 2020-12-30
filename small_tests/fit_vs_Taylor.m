%Simple test of Polynomial fitting versus Taylor expansion
%Recreate plots from https://math.libretexts.org/Bookshelves/Calculus/Book%3A_Active_Calculus_(Boelkins_et_al)/8%3A_Sequences_and_Series/8.5%3A_Taylor_Polynomials_and_Taylor_Series
close all;
clearvars;

N = 1000;               %Number of sample points of the sine
x = linspace(-5,5,N);   
y = sin(x); 
W = 0.25;               %Fit window (fit between -W<x<W)

%% Calculate 1st up to 9th order Taylor expansions of the sine
P1 = x;
P3 = x - x.^3/factorial(3);
P5 = x - x.^3/factorial(3) + x.^5/factorial(5);
P7 = x - x.^3/factorial(3) + x.^5/factorial(5) - x.^7/factorial(7);
P9 = x - x.^3/factorial(3) + x.^5/factorial(5) - x.^7/factorial(7) + x.^9/factorial(9);

%% Calculate 1st up to 9th order fit in windows W
%Determine which x/y values are eligible for fitting
xin = x(abs(x)<W);
yin = y(abs(x)<W);
%Determine coefficitions using polyfit
p1 = polyfit(xin,yin,1);
p3 = polyfit(xin,yin,3);
p5 = polyfit(xin,yin,5);
p7 = polyfit(xin,yin,7);
p9 = polyfit(xin,yin,9);
%Calculate values elsewhere using polyval
Q1 = polyval(p1,x);
Q3 = polyval(p3,x);
Q5 = polyval(p5,x);
Q7 = polyval(p7,x);
Q9 = polyval(p9,x);

%% Print maximum absolute error between Taylor and fit
fprintf('Max error fit - Taylor (1st order): %.2e\n',max(abs(P1-Q1)));
fprintf('Max error fit - Taylor (3rd order): %.2e\n',max(abs(P3-Q3)));
fprintf('Max error fit - Taylor (5th order): %.2e\n',max(abs(P5-Q5)));
fprintf('Max error fit - Taylor (7th order): %.2e\n',max(abs(P7-Q7)));
fprintf('Max error fit - Taylor (9th order): %.2e\n',max(abs(P9-Q9)));

%% Plot the results for 1st - 9th order Taylor
figure(1);
subplot(2,3,1);
plot(x,y),hold on;
plot(x,P1);
title('1st order Taylor');
ylim([-3 3]);
axis square;
grid on;

subplot(2,3,2);
plot(x,y),hold on;
plot(x,P3);
title('3rd order Taylor');
ylim([-3 3]);
axis square;
grid on;

subplot(2,3,3);
plot(x,y),hold on;
plot(x,P5);
title('5th order Taylor');
ylim([-3 3]);
axis square;
grid on;

subplot(2,3,4);
plot(x,y),hold on;
plot(x,P7);
title('7th order Taylor');
ylim([-3 3]);
axis square;
grid on;

subplot(2,3,5)
plot(x,y),hold on;
plot(x,P9);
title('9th order Taylor');
ylim([-3 3]);
axis square;
grid on;

%% Plot the results for 1st - 9th order polynomial fit
figure(2);
subplot(2,3,1);
plot(x,y),hold on;
plot(x,Q1);
title('1st order Fit');
ylim([-3 3]);
axis square;
grid on;

subplot(2,3,2);
plot(x,y),hold on;
plot(x,Q3);
title('3rd order Fit');
ylim([-3 3]);
axis square;
grid on;

subplot(2,3,3);
plot(x,y),hold on;
plot(x,Q5);
title('5th order Fit');
ylim([-3 3]);
axis square;
grid on;

subplot(2,3,4);
plot(x,y),hold on;
plot(x,Q7);
title('7th order Fit');
ylim([-3 3]);
axis square;
grid on;

subplot(2,3,5)
plot(x,y),hold on;
plot(x,Q9);
title('9th order Fit');
ylim([-3 3]);
axis square;
grid on;