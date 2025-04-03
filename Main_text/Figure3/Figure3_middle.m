%% ID oscillators

%well
format long

%mr clean
clc

%starting spot
M = 0;

%choose diagonal row
N = 1;

%cutoff point to approximate infinite continued fraction as a finite one
cutoff = 200;

%fixed system parameters 
D = 0.1;
w = 2;
tau = 0;

%coupling strength
gamma = 0:.0001:.02;

%solution vectors
lambda_1 = zeros(1,length(gamma));
lambda_2 = zeros(1,length(gamma));

%loop over all parameter space
for ii = 1:length(gamma)
        
    %initialize matrix
    my_matrix_friend = zeros(2*cutoff+1,2*cutoff+1);
    
    %form the matrix
    count = 1;
    for j=-cutoff+M:cutoff+M
        
        %diagonal entries
        my_matrix_friend(count,count) = Q_mid(N,j,D,w+tau,w);
        
        %super/subdiagonal entries
        my_matrix_friend(count,count+1) = Q_plus(N,j,gamma(ii));
        my_matrix_friend(count+1,count) = Q_minus(N,j+1,gamma(ii));
        
        %update count
        count = count+1;
    end
    
    if mod(N,2) == 0
        my_matrix_friend = my_matrix_friend(1:2*cutoff+1,1:2*cutoff+1);
    else
        my_matrix_friend(end,end) = Q_mid(N,cutoff+M+1,D,w+tau,w);
    end
    
    %diagonalize
    lambda_boi = eig(my_matrix_friend);
    
    %sort
    [~,temp] = sort(-real(lambda_boi));
    lambda_boi = lambda_boi(temp);
    
    %get the eigenvalues to compare
    lambda_1(ii) = lambda_boi(1);
    lambda_2(ii) = lambda_boi(2);
    if imag(lambda_1(ii)) > imag(lambda_2(ii))
        temp = lambda_2(ii);
        lambda_2(ii) = lambda_1(ii);
        lambda_1(ii) = temp;
    end         
end

%exact expressions
l1 = (-D+sqrt(-1)*w) + tau/2*sqrt(-1)  + 1/2*sqrt(gamma.^2-tau^2);
l1n = (-D+sqrt(-1)*w) + tau/2*sqrt(-1) - 1/2*sqrt(gamma.^2-tau^2);

%im part
figure(1)
subplot(2,1,1)
hold on
plot(gamma,imag(lambda_1),'.','color','k','markersize',30)
plot(gamma,imag(lambda_2),'.','color','k','markersize',30)
plot(gamma,imag(l1),'.','color',[.6 .6 .6 .6],'markersize',10)
plot(gamma,imag(l1n),'.','color',[.6 .6 .6 .6],'markersize',10)
ylabel('Im(\lambda)')
box on
set(gca,'fontsize',15)
%xline(0.335961852,'m--','linewidth',2)
%xline(0,'--','linewidth',2)
ylim([1.99 2.01])

subplot(2,1,2)
hold on
plot(gamma,real(lambda_1),'.','color','k','markersize',30)
plot(gamma,real(lambda_2),'.','color','k','markersize',30)
plot(gamma,real(l1),'.','color',[.6 .6 .6 .6],'markersize',10)
plot(gamma,real(l1n),'.','color',[.6 .6 .6 .6],'markersize',10)
xlabel('\kappa')
ylabel('Re(\lambda)')
box on
set(gca,'fontsize',15)
set(gcf,'position',[406,115,377,641])
%xline(0.335961852,'m--','linewidth',2)
%xline(0,'--','linewidth',2)


%% nonID oscillators

%starting spot
M = 0;

%choose diagonal row
N = 1;

%cutoff point to approximate infinite continued fraction as a finite one
cutoff = 200;

%fixed system parameters 
D = 0.1;
w = 2;
tau = .01;

%coupling strength
gamma = 0:.0001:.02;

%solution vectors
lambda_1 = zeros(1,length(gamma));
lambda_2 = zeros(1,length(gamma));

%loop over all parameter space
for ii = 1:length(gamma)
        
    %initialize matrix
    my_matrix_friend = zeros(2*cutoff+1,2*cutoff+1);
    
    %form the matrix
    count = 1;
    for j=-cutoff+M:cutoff+M
        
        %diagonal entries
        my_matrix_friend(count,count) = Q_mid(N,j,D,w+tau,w);
        
        %super/subdiagonal entries
        my_matrix_friend(count,count+1) = Q_plus(N,j,gamma(ii));
        my_matrix_friend(count+1,count) = Q_minus(N,j+1,gamma(ii));
        
        %update count
        count = count+1;
    end
    
    if mod(N,2) == 0
        my_matrix_friend = my_matrix_friend(1:2*cutoff+1,1:2*cutoff+1);
    else
        my_matrix_friend(end,end) = Q_mid(N,cutoff+M+1,D,w+tau,w);
    end
    
    %diagonalize
    lambda_boi = eig(my_matrix_friend);
    
    %sort
    [~,temp] = sort(-real(lambda_boi));
    lambda_boi = lambda_boi(temp);
    
    %get the eigenvalues to compare
    lambda_1(ii) = lambda_boi(1);
    lambda_2(ii) = lambda_boi(2);
    if imag(lambda_1(ii)) > imag(lambda_2(ii))
        temp = lambda_2(ii);
        lambda_2(ii) = lambda_1(ii);
        lambda_1(ii) = temp;
    end         
end

%exact expressions
l1 = (-D+sqrt(-1)*w) + tau/2*sqrt(-1)  + 1/2*sqrt(gamma.^2-tau^2);
l1n = (-D+sqrt(-1)*w) + tau/2*sqrt(-1) - 1/2*sqrt(gamma.^2-tau^2);

%im part
figure(2)
subplot(2,1,1)
hold on
plot(gamma,imag(lambda_1),'.','color','k','markersize',30)
plot(gamma,imag(lambda_2),'.','color','k','markersize',30)
plot(gamma,imag(l1),'.','color',[.6 .6 .6 .6],'markersize',10)
plot(gamma,imag(l1n),'.','color',[.6 .6 .6 .6],'markersize',10)
xline(0.01,'--','linewidth',2)
ylabel('Im(\lambda)')
box on
set(gca,'fontsize',15)
%xline(0.335961852,'m--','linewidth',2)
ylim([1.99 2.02])

subplot(2,1,2)
hold on
plot(gamma,real(lambda_1),'.','color','k','markersize',30)
plot(gamma,real(lambda_2),'.','color','k','markersize',30)
plot(gamma,real(l1),'.','color',[.6 .6 .6 .6],'markersize',10)
plot(gamma,real(l1n),'.','color',[.6 .6 .6 .6],'markersize',10)
xline(0.01,'--','linewidth',2)
xlabel('\kappa')
ylabel('Re(\lambda)')
box on
set(gca,'fontsize',15)
set(gcf,'position',[406,115,377,641])
%xline(0.335961852,'m--','linewidth',2)


%% functionals to make my life easier

%Q-
function[u] = Q_minus(N,m,gamma)
u = gamma/2*(N-2*(m-1));
end

%Q
function[u] = Q_mid(N,m,D,w1,w2)
u = -D*(m^2+(N-m)^2)+sqrt(-1)*w1*m+sqrt(-1)*w2*(N-m);
end

%Q+
function[u] = Q_plus(N,m,gamma)
u = gamma/2*(2*(m+1)-N);
end
