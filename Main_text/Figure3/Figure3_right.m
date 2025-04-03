%% ID osc

%mr clean
clc
clf

%parameters
w = 1;
tau = 0;

%coupling
gamma = 0:.0001:.02; 

%size of mah boi
n = 3;

%solution vector of eigenvalues
LAMBDA = zeros(2,length(gamma));

%loopin' it
for pp = 1:length(gamma)
    
    %joint matrix without coupling
    A = -eye(n)*(w+tau) + (w+tau)*circshift(eye(n),-1);
    B = -eye(n)*w + w*circshift(eye(n),-1);
    well = kron(B,eye(n)) + kron(eye(n),A);

    %well
    GAMMA = gamma(pp);
    
    %form the coupling structure vector
    if mod(n,2) == 0
        g = [0, GAMMA*ones(1,(n/2)-1), 0, -GAMMA*ones(1,(n/2)-1)];
    else
        g = [0, GAMMA*ones(1,((n-1)/2)), -GAMMA*ones(1,((n-1)/2))];
    end
    
    %form coupling (for beta)
    I = eye(n);
    I = kron(circshift(I,-1),I);
    G = [];
    for i=0:n-1
        G = [G circshift(g,i)];
    end
    BETA = I.*G';
    
    %form coupling (for alpha)
    I = eye(n);
    I = kron(I,circshift(I,-1));
    G = [];
    for i=0:n-1
        G = [G circshift(-g,i)];
    end
    ALPHA = I.*G';
       
    % Update the joint transition matrix with coupling
    well = (well + BETA + ALPHA)*1;
    
    %diagonalize
    [~,lambda] = eig(well);
    lambda = diag(lambda);
    lambda = lambda(imag(lambda)>0);
    [~,I] = sort(real(lambda),'descend');
    lambda = lambda(I);
    LAMBDA(:,pp) = lambda(1:2);
end

%exact solutions (to check)
l1 = (-3/2+sqrt(3)/2*sqrt(-1))*w + (3/4 - sqrt(3)/4*sqrt(-1))*(-tau) + sqrt(2)/4*sqrt(-gamma.^2*(4-2*sqrt(12)*sqrt(-1))+tau^2*(3-3*sqrt(12)/2*sqrt(-1)));
l2 = (-3/2+sqrt(3)/2*sqrt(-1))*w + (3/4 - sqrt(3)/4*sqrt(-1))*(-tau) - sqrt(2)/4*sqrt(-gamma.^2*(4-2*sqrt(12)*sqrt(-1))+tau^2*(3-3*sqrt(12)/2*sqrt(-1)));

%im part
figure(1)
subplot(2,1,1)
hold on
plot(gamma,imag(LAMBDA(1,:)),'.','color','k','markersize',30)
plot(gamma,imag(LAMBDA(2,:)),'.','color','k','markersize',30)
plot(gamma,imag(l1),'.','color',[255, 182, 193 256]/256,'markersize',10);
plot(gamma,imag(l2),'.','color',[255, 182, 193 256]/256,'markersize',10);
ylabel('Im(\lambda)')
box on
set(gca,'fontsize',15)
%xline(0.335961852,'m--','linewidth',2)
%xline(0,'--','linewidth',2)
ylim([0.84 0.89])
%xlim([-.1 .1])

subplot(2,1,2)
hold on
plot(gamma,real(LAMBDA(1,:)),'.','color','k','markersize',30)
plot(gamma,real(LAMBDA(2,:)),'.','color','k','markersize',30)
plot(gamma,real(l1),'.','color',[255, 182, 193 256]/256,'markersize',10);
plot(gamma,real(l2),'.','color',[255, 182, 193 256]/256,'markersize',10);
xlabel('\kappa')
ylabel('Re(\lambda)')
box on
set(gca,'fontsize',15)
set(gcf,'position',[406,115,377,641])
%xline(0.335961852,'m--','linewidth',2)
%xline(0,'--','linewidth',2)
%xlim([-.1 .1])
ylim([-1.515 -1.485])


%% nonID osc

%parameters
w = 1;
tau = 0.01;

%coupling
gamma = 0:.0001:.02; 

%size of mah boi
n = 3;

%solution vector of eigenvalues
LAMBDA = zeros(2,length(gamma));

%loopin' it
for pp = 1:length(gamma)
    
    %joint matrix without coupling
    A = -eye(n)*(w+tau) + (w+tau)*circshift(eye(n),-1);
    B = -eye(n)*w + w*circshift(eye(n),-1);
    well = kron(B,eye(n)) + kron(eye(n),A);

    %well
    GAMMA = gamma(pp);
    
    %form the coupling structure vector
    if mod(n,2) == 0
        g = [0, GAMMA*ones(1,(n/2)-1), 0, -GAMMA*ones(1,(n/2)-1)];
    else
        g = [0, GAMMA*ones(1,((n-1)/2)), -GAMMA*ones(1,((n-1)/2))];
    end
    
    %form coupling (for beta)
    I = eye(n);
    I = kron(circshift(I,-1),I);
    G = [];
    for i=0:n-1
        G = [G circshift(g,i)];
    end
    BETA = I.*G';
    
    %form coupling (for alpha)
    I = eye(n);
    I = kron(I,circshift(I,-1));
    G = [];
    for i=0:n-1
        G = [G circshift(-g,i)];
    end
    ALPHA = I.*G';
       
    % Update the joint transition matrix with coupling
    well = (well + BETA + ALPHA)*1;
    
    %diagonalize
    [~,lambda] = eig(well);
    lambda = diag(lambda);
    lambda = lambda(imag(lambda)>0);
    [~,I] = sort(real(lambda),'descend');
    lambda = lambda(I);
    LAMBDA(:,pp) = lambda(1:2);
end

%exact solutions (to check)
l1 = (-3/2+sqrt(3)/2*sqrt(-1))*w + (3/4 - sqrt(3)/4*sqrt(-1))*(-tau) + sqrt(2)/4*sqrt(-gamma.^2*(4-2*sqrt(12)*sqrt(-1))+tau^2*(3-3*sqrt(12)/2*sqrt(-1)));
l2 = (-3/2+sqrt(3)/2*sqrt(-1))*w + (3/4 - sqrt(3)/4*sqrt(-1))*(-tau) - sqrt(2)/4*sqrt(-gamma.^2*(4-2*sqrt(12)*sqrt(-1))+tau^2*(3-3*sqrt(12)/2*sqrt(-1)));

%im part
figure(2)
subplot(2,1,1)
hold on
plot(gamma,imag(LAMBDA(1,:)),'.','color','k','markersize',30)
plot(gamma,imag(LAMBDA(2,:)),'.','color','k','markersize',30)
plot(gamma,imag(l1),'.','color',[255, 182, 193 256]/256,'markersize',10);
plot(gamma,imag(l2),'.','color',[255, 182, 193 256]/256,'markersize',10);
ylabel('Im(\lambda)')
box on
set(gca,'fontsize',15)
%xline(0.335961852,'m--','linewidth',2)
xline(sqrt(3)/2*tau,'--','linewidth',2)
%ylim([1.9 2.1])
%xlim([-.1 .1])

subplot(2,1,2)
hold on
plot(gamma,real(LAMBDA(1,:)),'.','color','k','markersize',30)
plot(gamma,real(LAMBDA(2,:)),'.','color','k','markersize',30)
plot(gamma,real(l1),'.','color',[255, 182, 193 256]/256,'markersize',10);
plot(gamma,real(l2),'.','color',[255, 182, 193 256]/256,'markersize',10);
xlabel('\kappa')
ylabel('Re(\lambda)')
box on
set(gca,'fontsize',15)
set(gcf,'position',[406,115,377,641])
%xline(0.335961852,'m--','linewidth',2)
xline(sqrt(3)/2*tau,'--','linewidth',2)
%xlim([-.1 .1])










