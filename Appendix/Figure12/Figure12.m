%% setup

%mr clean
clc

%parameters
a = 1*sqrt(3);
b = 2*sqrt(3);

%coupling
gamma_star = sqrt(3)/2*abs(b-a);
gamma = 0:.01:1.7;

%size
n = 3;

%phase and abs vecs
Absa = zeros(3,length(gamma));
Absb = Absa;
ABSa = zeros(3,length(gamma));
ABSb = ABSa;

%full phase difference
phase_dif = zeros(9,length(gamma));

%quality factor
QF = zeros(2,length(gamma));

%loopin
for pp=1:length(gamma)
    
    %% form the matrix
    
    %joint matrix without coupling
    A = -eye(n)*a + a*circshift(eye(n),-1);
    B = -eye(n)*b + b*circshift(eye(n),-1);
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
    well = (well + BETA + ALPHA).';
    
    
    %% diagonalize
    
    %diagonalize
    [w,lambda] = eig(well.');
    lambda = diag(lambda);
    big_lambda = lambda;
    
    %sort
    index = imag(lambda)>0;
    w = w(:,index);
    lambda = lambda(index);
    
    [~,I] = sort(real(lambda),'descend');
    lambda = lambda(I);
    w = w(:,I);
    
    %form the eigenvalues and eigenvectors for each oscillator
    lambda1 = lambda(1);
    lambda2 = lambda(2);
    
    w1 = w(:,1);
    w2 = w(:,2);
    
    %record quality factor
    QF(1,pp) = abs(imag(lambda1))/abs(real(lambda1));
    QF(2,pp) = abs(imag(lambda2))/abs(real(lambda2));
    
    
    %% projection matrices
    
    %form them
    Pa = repmat(eye(n),[1, n]);
    Pb = kron(eye(n),ones(1,n));
    
    
    %% project the heck out of those rascals
    
    %project
    PQa = Pa*w1;
    PQb = Pb*w2;
    
    %normalize
    [~,index_a] = (max(abs(real(PQa))));
    PQa = PQa/-sign(real(PQa(index_a)));
    [~,index_b] = (max(abs(real(PQb))));
    PQb = PQb/-sign(real(PQb(index_b)));
    
    %test
    Absa(:,pp) = angle(PQa)+pi;
    Absb(:,pp) = angle(PQb)+pi;
    ABSa(:,pp) = abs(PQa);
    ABSb(:,pp) = abs(PQb);
    
    %full phase difference
    phase_dif(:,pp) = angle(w2) - angle(w1); 
end

%% welp

%do stuff to the 1st oscillator phase for data analysis
Absa(Absa > 4) = Absa(Absa > 4) - 2*pi/3;
Absa(Absa > 2) = Absa(Absa > 2) - 2*pi/3;
Absa(Absa > 1) = Absa(Absa > 1) - pi/3;
Absa(2,:) = Absa(1,:) + 2*pi/3;
Absa(3,:) = Absa(2,:) + 2*pi/3;

%do stuff to the 1st oscillator phase for data analysis
Absb(Absb<1) = Absb(Absb<1) + 2*pi;
Absb(Absb>5) = Absb(Absb>5) - 2*pi/3;
Absb(Absb>3) = Absb(Absb>3) - 2*pi/3;
Absb = Absb - 2*pi/3;
Absb(2,:) = Absb(1,:) + 2*pi/3;
Absb(3,:) = Absb(2,:) + 2*pi/3;

figure(5)
hold on
plot(gamma,Absa(1,:)-Absb(1,:),'k.','markersize',40)
plot(gamma,Absa(2,:)-Absb(2,:),'k.','markersize',30)
plot(gamma,Absa(3,:)-Absb(3,:),'k.','markersize',20)
xline(gamma_star,'--','linewidth',2)
%legend('Component 1','Component 2','Component 3','bifurcation pt')
xlabel('\kappa')
ylabel('Arg(P_AQ_1^*) - Arg(P_BQ_2^*)')
set(gca,'fontsize',12)
box on
axis square
grid on
