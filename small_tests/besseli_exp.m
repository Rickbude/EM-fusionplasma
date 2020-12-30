clearvars;
close all;

W          = 2;
rho_L      = 1e-6;
k_perp_max = W/rho_L;
k_perp     = linspace(0,k_perp_max,1000);

lambda = 0.5*k_perp.^2.*rho_L.^2;

Lambda_0 = besseli(0,lambda).*exp(-lambda);
Lambda_1 = besseli(1,lambda).*exp(-lambda);
Lambda_2 = besseli(2,lambda).*exp(-lambda);

figure(1);
plot(lambda,Lambda_0), hold on
plot(lambda,Lambda_1), hold on
plot(lambda,Lambda_2), hold on

xlabel('\lambda')
ylabel('\Lambda_n')
legend('\Lambda_0','\Lambda_1','\Lambda_2');

figure(2);
subplot(1,2,1);
plot(lambda,Lambda_0,'Linewidth',2),hold on;

approx_0_4 = lambda.^0*d_Lambda_n(0,0,0)/factorial(0) + ...
           lambda.^1*d_Lambda_n(1,0,0)/factorial(1) + ...
           lambda.^2*d_Lambda_n(2,0,0)/factorial(2) + ...
           lambda.^3*d_Lambda_n(3,0,0)/factorial(3) + ...
           lambda.^4*d_Lambda_n(4,0,0)/factorial(4);       
       
fit_0_4    = polyfit(lambda,Lambda_0,4);

plot(lambda,approx_0_4, '-.','Linewidth',2);
plot(lambda,polyval(fit_0_4,lambda),':','Linewidth',2)

xlabel('\lambda')
ylabel('\Lambda_0')
legend('\Lambda_0','Taylor (4th order)','fit (4th order)');

subplot(1,2,2);
approx_0_8 = lambda.^0*d_Lambda_n(0,0,0)/factorial(0) + ...
             lambda.^1*d_Lambda_n(1,0,0)/factorial(1) + ...
             lambda.^2*d_Lambda_n(2,0,0)/factorial(2) + ...
             lambda.^3*d_Lambda_n(3,0,0)/factorial(3) + ...
             lambda.^4*d_Lambda_n(4,0,0)/factorial(4) + ...
             lambda.^5*d_Lambda_n(5,0,0)/factorial(5) + ...
             lambda.^6*d_Lambda_n(6,0,0)/factorial(6) + ...
             lambda.^7*d_Lambda_n(7,0,0)/factorial(7) + ...
             lambda.^8*d_Lambda_n(8,0,0)/factorial(8);       
       
fit_0_8    = polyfit(lambda,Lambda_0,8);


plot(lambda,Lambda_0,'Linewidth',2),hold on;
plot(lambda,approx_0_8, '-.','Linewidth',2);
plot(lambda,polyval(fit_0_8,lambda),':','Linewidth',2)

xlabel('\lambda')
ylabel('\Lambda_0')
legend('\Lambda_0','Taylor (8th order)','fit (8th order)');


figure(4);

approx_1_10= lambda.^0*d_Lambda_n(0,1,0)/factorial(0) + ...
             lambda.^1*d_Lambda_n(1,1,0)/factorial(1) + ...
             lambda.^2*d_Lambda_n(2,1,0)/factorial(2) + ...
             lambda.^3*d_Lambda_n(3,1,0)/factorial(3) + ...
             lambda.^4*d_Lambda_n(4,1,0)/factorial(4) + ...
             lambda.^5*d_Lambda_n(5,1,0)/factorial(5) + ...
             lambda.^6*d_Lambda_n(6,1,0)/factorial(6) + ...
             lambda.^7*d_Lambda_n(7,1,0)/factorial(7) + ...
             lambda.^8*d_Lambda_n(8,1,0)/factorial(8) + ...
             lambda.^9*d_Lambda_n(9,1,0)/factorial(9) + ...
             lambda.^10*d_Lambda_n(10,1,0)/factorial(10);       
       
fit_1_10    = polyfit(lambda,Lambda_1,10);


plot(lambda,Lambda_1,'Linewidth',2),hold on;
plot(lambda,approx_1_10, '-.','Linewidth',2);
plot(lambda,polyval(fit_1_10,lambda),':','Linewidth',2)

xlabel('\lambda')
ylabel('\Lambda_1')
legend('\Lambda_1','Taylor (10th order)','fit (10th order)');

figure(5)
plot(lambda,Lambda_0,'Linewidth',2),hold on;

approx_0_4 = lambda.^0*d_Lambda_n(0,0,0)/factorial(0) + ...
           lambda.^1*d_Lambda_n(1,0,0)/factorial(1) + ...
           lambda.^2*d_Lambda_n(2,0,0)/factorial(2) + ...
           lambda.^3*d_Lambda_n(3,0,0)/factorial(3) + ...
           lambda.^4*d_Lambda_n(4,0,0)/factorial(4);       
       
fit_0_4    = polyfit(lambda,Lambda_0,4);

plot(lambda,approx_0_4, '-.','Linewidth',2);
plot(lambda,polyval(fit_0_4,lambda),':','Linewidth',2)

xlabel('\lambda')
ylabel('I_0(\lambda)e^{-\lambda}')
legend('I_0(\lambda)e^{-\lambda}','Taylor (4th order)','fit (4th order)','location','northwest');
