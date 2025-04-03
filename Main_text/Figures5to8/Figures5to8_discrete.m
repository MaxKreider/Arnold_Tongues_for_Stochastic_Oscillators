%% setup

%mr clean, gravy why you flow so mean
clc
clf

%system parameters
alpha = 10;
beta = 15;
kappa = 5*sqrt(3)/2;
kappa = 6;

%number of trials
M = 1000;

%interpolation time vector and Tmax
Delta = 1/10000;
N = 2^15;
t_query = 0:Delta:(N-1)*Delta;
Tmax = t_query(end) + 10000*Delta;

%frequency vector
kkk = (-N/2:N/2-1);
fk = 1/(N*Delta)*kkk*2*pi;

%power spectra
Plord_a = zeros(1,length(fk));
Plord_b = zeros(1,length(fk));
Plord1 = zeros(1,length(fk));
Plord2 = zeros(1,length(fk));
CrossLord = zeros(1,length(fk));
CrossLord_og = zeros(1,length(fk));


%% construct the state transition matrix

%columns
a1 = [-alpha-beta beta 0 alpha 0 0 0 0 0]';
a2 = [0 -alpha-beta beta-kappa 0 alpha+kappa 0 0 0 0]';
a3 = [beta+kappa 0 -alpha-beta 0 0 alpha-kappa 0 0 0]';
a4 = [0 0 0 -alpha-beta beta+kappa 0 alpha-kappa 0 0]';
a5 = [0 0 0 0 -alpha-beta beta 0 alpha 0]';
a6 = [0 0 0 beta-kappa 0 -alpha-beta 0 0 alpha+kappa]';
a7 = [alpha+kappa 0 0 0 0 0 -alpha-beta beta-kappa 0]';
a8 = [0 alpha-kappa 0 0 0 0 0 -alpha-beta beta+kappa]';
a9 = [0 0 alpha 0 0 0 beta 0 -alpha-beta]';

%form the matrix
A = [a1 a2 a3 a4 a5 a6 a7 a8 a9];


%% find the Q functions

%Q functions and eigenvalues
[w,lambda] = eig(A.');
lambda = diag(lambda);

%sort the eigenvalues by keeping only those with positive imaginary parts
I = imag(lambda)>0;
lambda = lambda(I);
w = w(:,I);

%further sort the eigenvalues by real parts
[~,I] = sort(real(lambda),'descend');
lambda = lambda(I);
w = w(:,I);

%output to debug
lambda1 = lambda(1)
lambda2 = lambda(2)
Q1 = conj(w(:,1))
Q2 = conj(w(:,2))


%% find the stationary distribution

%find it
[v,lambda] = eig(A);
lambda = diag(lambda);

%sort according to real parts
[~,I] = sort(real(lambda),'descend');
lambda = lambda(I);
v = v(:,I);
P0 = v(:,1)/sum(v(:,1));

%for choosing a random initial condition
P0_vec = cumsum(P0);


%% form the stoichiometry matrix

%columns
s1 = [-1 1 0 0 0 0 0 0 0]';
s2 = [-1 0 0 1 0 0 0 0 0]';
s3 = [0 -1 1 0 0 0 0 0 0]';
s4 = [0 -1 0 0 1 0 0 0 0]';
s5 = [1 0 -1 0 0 0 0 0 0]';
s6 = [0 0 -1 0 0 1 0 0 0]';
s7 = [0 0 0 -1 1 0 0 0 0]';
s8 = [0 0 0 -1 0 0 1 0 0]';
s9 = [0 0 0 0 -1 1 0 0 0]';
s10 = [0 0 0 0 -1 0 0 1 0]';
s11 = [0 0 0 1 0 -1 0 0 0]';
s12 = [0 0 0 0 0 -1 0 0 1]';
s13 = [0 0 0 0 0 0 -1 1 0]';
s14 = [1 0 0 0 0 0 -1 0 0]';
s15 = [0 0 0 0 0 0 0 -1 1]';
s16 = [0 1 0 0 0 0 0 -1 0]';
s17 = [0 0 0 0 0 0 1 0 -1]';
s18 = [0 0 1 0 0 0 0 0 -1]';

%form this rascal
S = [s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 s14 s15 s16 s17 s18];


%% perform the gillespie algorithm

%multiple trials
for mmm = 1:M
    
    %establish solution vectors for each run
    SOL = zeros(1,length(t_query));
    SOL_Q1 = zeros(1,length(t_query));
    SOL_Q2 = zeros(1,length(t_query));
    power_a = 0;
    power_b = 0;
    power_x_Q1 = 0;
    power_x_Q2 = 0;
    power_cross = 0;
    power_cross_og = 0;
    
    
    %% initial amount of each species (chosen at random)    
    
    %initialize IC vector
    x0 = zeros(9,1);
    
    %choose where to start
    xi = rand;
    count = 1;
    while xi > P0_vec(count)
        count = count + 1;
    end
    x0(count) = 9;
    
    
    %% define propensity function
    prop = @(x) [beta*x(1); alpha*x(1);
                (beta-kappa)*x(2); (alpha+kappa)*x(2);
                (beta+kappa)*x(3); (alpha-kappa)*x(3);
                (beta+kappa)*x(4); (alpha-kappa)*x(4);
                beta*x(5); alpha*x(5);
                (beta-kappa)*x(6); (alpha+kappa)*x(6);
                (beta-kappa)*x(7); (alpha+kappa)*x(7);
                (beta+kappa)*x(8); (alpha-kappa)*x(8);
                beta*x(9); alpha*x(9)];
    
    
    %% run the simulation
    
    %starting time of the simulation
    t=0;
    
    %initialize
    tspan=[];
    tspan(1)=t;
    X=[];
    X(:,1)=x0;
    
    %count
    count = 1;
    
    %solve
    while t<Tmax
        
        %draw the next reaction time
        bigalpha=prop(X(:,count));
        xi=rand;
        s=-1/sum(bigalpha)*log(1-xi);
        
        %draw the reaction
        phi=0;
        ell=0;
        eta=rand;
        while (phi < eta)
            ell=ell+1;
            phi=phi + bigalpha(ell)/sum(bigalpha);
        end
        
        %update
        t=t+s;
        X(:,count+1)=X(:,count) + S(:,ell);
        tspan(count+1)=t;
        
        %update count
        count = count + 1;
    end
    
    
    %% probably not a strictly correct thing to do, but here we are ;)
    
    %interpolate the data onto an evenly spaced time grid
    X_interped = interp1(tspan,X',t_query,'previous')';
        
    
    %% transform to Q function coordinates
    
    %loop over all time points
    for pp=1:length(t_query)
        
       %coordinate transform 
       SOL_Q1(pp) = dot(Q1,X_interped(:,pp));
       SOL_Q2(pp) = dot(Q2,X_interped(:,pp));
    end
    
    
    %% take the fourier transform and compute power/cross spectra
    
    %fft
    Fc_Q1 = fftshift(fft(SOL_Q1 - 0*mean(SOL_Q1)));
    Fc_Q2 = fftshift(fft(SOL_Q2 - 0*mean(SOL_Q2)));
    Fc_a = fftshift(fft(X_interped(1,:) - 0*mean(X_interped(1,:))));
    Fc_b = fftshift(fft(X_interped(8,:) - 0*mean(X_interped(8,:))));
        
    %power spectra
    power_a = power_a + (Fc_a).*conj(Fc_a);
    power_b = power_b + (Fc_b).*conj(Fc_b);
    power_cross_og = power_cross_og + conj(Fc_a).*(Fc_b);
    power_x_Q1 = power_x_Q1 + (Fc_Q1).*conj(Fc_Q1);
    power_x_Q2 = power_x_Q2 + (Fc_Q2).*conj(Fc_Q2);
    power_cross = power_cross + conj(Fc_Q1).*(Fc_Q2);
    
    
    %aggregate
    Plord_a = Plord_a + power_a;
    Plord_b = Plord_b + power_b;
    Plord1 = Plord1 + power_x_Q1;
    Plord2 = Plord2 + power_x_Q2;
    CrossLord = CrossLord + power_cross;
    CrossLord_og = CrossLord_og + power_cross_og;
end


%% numerical solutions

%normalize numerical solutions
Plord_a = Plord_a/(M*N)*Delta;
Plord_b = Plord_b/(M*N)*Delta;
Plord1 = Plord1/(M*N)*Delta;
Plord2 = Plord2/(M*N)*Delta;
CrossLord = CrossLord/(M*N)*Delta;
CrossLord_og = CrossLord_og/(M*N)*Delta;

%well
I=find(Plord_a>.4);
Plord_a(I) = (Plord_a(I-1)+Plord_a(I+1))/2;
I=find(Plord_b>.4);
Plord_b(I) = (Plord_b(I-1)+Plord_b(I+1))/2;
I=find(CrossLord_og>.4);
CrossLord_og(I) = (CrossLord_og(I-1)+CrossLord_og(I+1))/2;


%% exact solutions

%exact solutions for power spectra
exact_power_1 = -2*real(lambda1)./(real(lambda1)^2 + (fk-imag(lambda1)).^2);
exact_power_2 = -2*real(lambda2)./(real(lambda2)^2 + (fk-imag(lambda2)).^2);
exact_cross = -(Q2'*Q1)*(1./((lambda2)-sqrt(-1)*fk) + 1./(conj(lambda1)+sqrt(-1)*fk));


%% visualize

%plot power spectra
figure(5)
hold on
plot(fk,real(Plord_a),'-','color',[0.8500 0.3250 0.0980 .3],'linewidth',7)
plot(fk,real(Plord_b),'-','color',[0.4940 0.1840 0.5560 .3],'linewidth',7)
plot(fk,Plord1,'-','color',[0.8500 0.3250 0.0980 .9],'linewidth',10)
plot(fk,exact_power_1,'-','color','y','linewidth',3)
plot(fk,Plord2,'-','color',[0.4940 0.1840 0.5560 .9],'linewidth',10)
plot(fk,exact_power_2,'-','color','#FFC0CB','linewidth',3)
xlim([-60 60])
ylim([0 .2])
xlabel('frequency \nu')
ylabel('S_1(\nu)')
box on
axis square
set(gca,'fontsize',15)

%plot cross spectra
figure(6)
hold on
plot(fk,real(CrossLord_og),'-','color',[.6 .6 .6 .6],'linewidth',10)
plot(fk,real(CrossLord),'-','color','#006400','linewidth',10)
plot(fk,real(exact_cross),'-','color',[0 1 0],'linewidth',3)
xlim([-60 60])
ylim([-.2 .2])
xlabel('frequency \nu')
ylabel('Re($S_{\lambda_x,\lambda_y}$)', 'Interpreter', 'latex');
box on
axis square
set(gca,'fontsize',15)

figure(7)
hold on
plot(fk,imag(CrossLord_og),'-','color',[.6 .6 .6 .6],'linewidth',10)
plot(fk,imag(CrossLord),'-','color',[0 0.4470 0.7410 .9],'linewidth',10)
plot(fk,imag(exact_cross),'-','color',[0 1 1],'linewidth',3)
xlim([-60 60])
ylim([-.2 .2])
xlabel('frequency \nu')
ylabel('Im($S_{\lambda_x,\lambda_y}$)', 'Interpreter', 'latex');
box on
axis square
set(gca,'fontsize',15)







