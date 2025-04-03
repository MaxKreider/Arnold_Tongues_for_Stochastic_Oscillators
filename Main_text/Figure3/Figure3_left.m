%% nonID oscillators

%mr clean
clc
format

%parameters
w = 2;
tau = 0.01;
k = 0.1;

%critical coupling strength
gamma_star = abs(tau)/2;

%gamma vector
Gamma = 0:.0001:.02;

for pp=1:length(Gamma)
    gamma = Gamma(pp);
    
    %do it
    A = [-k w+tau   0  0;
        -w-tau -k  0  0;
        0    0  -k w;
        0    0  -w -k];
    
    B = [-gamma, 0, gamma, 0;
        0, -gamma, 0, gamma;
        gamma, 0, -gamma, 0;
        0, gamma, 0, -gamma];
    
    A = A + B;
    
    [~,d] = eig(A);
    di(1:4,pp) = diag(d);
end

%from exact eigenvalue expression
l1 = (-k+sqrt(-1)*w) + tau/2*sqrt(-1) - Gamma + 1/2*sqrt(4*Gamma.^2-tau^2);
l1n = (-k+sqrt(-1)*w) + tau/2*sqrt(-1) - Gamma - 1/2*sqrt(4*Gamma.^2-tau^2);

figure(1)
subplot(2,1,2)
xline(gamma_star,'--','linewidth',2)
hold on
plot(Gamma,real(di),'k.','markersize',30)
plot(Gamma,real(l1),'g.','markersize',10)
plot(Gamma,real(l1n),'g.','markersize',10)
set(gca,'fontsize',15)
xlabel('\kappa')
box on
ylabel('Re(\lambda)')
ylim([-0.14 -0.09])

subplot(2,1,1)
hold on
xline(gamma_star,'--','linewidth',2)
plot(Gamma,imag(di),'k.','markersize',30)
plot(Gamma,imag(l1),'g.','markersize',10)
plot(Gamma,imag(l1n),'g.','markersize',10)
ylim([1.99 2.02])
set(gca,'fontsize',15)
box on
ylabel('Im(\lambda)')

set(gcf,'position',[406,115,377,641])


%% ID oscillators

%parameters
w = 2;
tau = 0;
k = 0.1;

%critical coupling strength
gamma_star = abs(tau)/2;

%gamma vector
Gamma = 0:.0001:.02;

for pp=1:length(Gamma)
    gamma = Gamma(pp);
    
    %do it
    A = [-k w+tau   0  0;
        -w-tau -k  0  0;
        0    0  -k w;
        0    0  -w -k];
    
    B = [-gamma, 0, gamma, 0;
        0, -gamma, 0, gamma;
        gamma, 0, -gamma, 0;
        0, gamma, 0, -gamma];
    
    A = A + B;
    
    [~,d] = eig(A);
    di(1:4,pp) = diag(d);
end

%from exact eigenvalue expression
l1 = (-k+sqrt(-1)*w) + tau/2*sqrt(-1) - Gamma + 1/2*sqrt(4*Gamma.^2-tau^2);
l1n = (-k+sqrt(-1)*w) + tau/2*sqrt(-1) - Gamma - 1/2*sqrt(4*Gamma.^2-tau^2);

figure(2)
subplot(2,1,2)
%xline(gamma_star,'--','linewidth',2)
hold on
plot(Gamma,real(di),'k.','markersize',30)
plot(Gamma,real(l1),'g.','markersize',10)
plot(Gamma,real(l1n),'g.','markersize',10)
set(gca,'fontsize',15)
xlabel('\kappa')
box on
ylabel('Re(\lambda)')

subplot(2,1,1)
hold on
%xline(gamma_star,'--','linewidth',2)
plot(Gamma,imag(di),'k.','markersize',30)
plot(Gamma,imag(l1),'g.','markersize',10)
plot(Gamma,imag(l1n),'g.','markersize',10)
ylim([1.99 2.01])
set(gca,'fontsize',15)
box on
ylabel('Im(\lambda)')

set(gcf,'position',[406,115,377,641])



