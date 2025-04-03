%% setup

%format
format long

%mr clean
clc

%system parameters
w1 = 2;
w2 = 2.5; %w2 = w + tau
D = 0.1;

%variable system parameters
kappa = 0.335961852;
kappa = 0.38;

%setup for psd do-da
Delta = 1/100;
N = 2^16;   
t = 0:Delta:(N-1)*Delta;

%frequency vector
k = (-N/2:N/2-1);
fk = 1/(N*Delta)*k*2*pi;

%number of trials
M = 20000;
M = 2500;


%% do it

%compute eigenvalues
LAMBDA = compute_eig(kappa,D,w1,w2);

%compute eigenfunction weights
[cfx,cbx] = compute_efunc(1,kappa,D,w1,w2,LAMBDA(1));
[cfy,cby] = compute_efunc(0,kappa,D,w1,w2,LAMBDA(2));

%compute Q function at that point
h = 0.01;
x = 0:h:2*pi;
y = x';
Q1 = computeQ(x,y,cfx,cbx,1);
Q2 = computeQ(x,y,cfy,cby,0);

%compute stationary distribution
[cfP0,cbP0] = compute_efunc_forward(0,kappa,D,w1,w2,0);
P0 = computeP0(x,y,cfP0,cbP0,0);
P0 = P0/trapz(y,trapz(x,P0,2));

%normalize Q functions to have unit variance
IQ1 = trapz(y,trapz(x,abs(Q1).^2.*P0,2));
IQ2 = trapz(y,trapz(x,abs(Q2).^2.*P0,2));
Q1 = Q1/sqrt(IQ1);
Q2 = Q2/sqrt(IQ2);

%rotate Q2 according to criteria
factor = (log(trapz(y,trapz(x,conj(Q2).*Q1.*P0,2))/trapz(y,trapz(x,conj(Q1).*Q2.*P0,2))));
if imag(factor) < 0
    factor = factor + 2*pi*sqrt(-1);
end
alpha_star = -sqrt(-1)/2*factor;
Q2 = Q2*exp(sqrt(-1)*alpha_star);

%computed bracketed mean term
I1 = trapz(y,trapz(x,Q1.*conj(Q2).*P0,2))

%yes
power_x = 0;
power_y = 0;
power_cross = 0;
dad_x = 0;
dad_y = 0;
dad_cross = 0;
    

%% compute numerical trajectories

parfor qq=1:M
        
    %solution vector
    Qx = zeros(1,length(t));
    Qy = zeros(1,length(t));
    XX = zeros(1,length(t));
    YY = zeros(1,length(t));
    
    %trajectory initial value (uniformly random on the torus)
    SOL = rand(2,1)*2*pi;
    
    %initial perturbed phase
    Qx(1) = computeQ(SOL(1),SOL(2),cfx,cbx,1);
    Qy(1) = computeQ(SOL(1),SOL(2),cfy,cby,0);
    XX(1) = cos(SOL(1)) + sqrt(-1)*sin(SOL(1));
    YY(1) = cos(SOL(2)) + sqrt(-1)*sin(SOL(2));
    
    %compute phases
    for pp=1:length(t)-1
        
        %EM method
        SOL = SOL + Delta*[w1+kappa*sin(SOL(2)-SOL(1));w2+kappa*sin(SOL(1)-SOL(2))] + sqrt(Delta)*sqrt(2*D)*randn(2,1);
        
        %compute Q function at that point
        Qx(pp+1) = computeQ(SOL(1),SOL(2),cfx,cbx,1)/sqrt(IQ1);
        Qy(pp+1) = computeQ(SOL(1),SOL(2),cfy,cby,0)/sqrt(IQ2)*exp(sqrt(-1)*alpha_star);
        XX(pp+1) = cos(SOL(1)) + sqrt(-1)*sin(SOL(1));
        YY(pp+1) = cos(SOL(2)) + sqrt(-1)*sin(SOL(2));
        
    end
    
    %take fft
    Fc_x = fftshift(fft(Qx));
    Fc_y = fftshift(fft(Qy));
    Fc_X = fftshift(fft(XX));
    Fc_Y = fftshift(fft(YY));
    
    %power spectra
    power_x = power_x + abs(Fc_x).^2;
    power_y = power_y + abs(Fc_y).^2;
    dad_x = dad_x + abs(Fc_X).^2;
    dad_y = dad_y + abs(Fc_Y).^2;
    
    %cross spectra
    power_cross = power_cross  + (Fc_x).*conj(Fc_y);
    dad_cross = dad_cross + (Fc_X).*conj(Fc_Y);
    
end

%normalize
power_x = power_x/M;
power_y = power_y/M;
power_cross = power_cross/M;
dad_x = dad_x/M;
dad_y = dad_y/M;
dad_cross = dad_cross/M;

%exact solutions for power spectra
mu_x = real(LAMBDA(1));
omega_x = imag(LAMBDA(1));
mu_y = real(LAMBDA(2));
omega_y = imag(LAMBDA(2));
exact_x = 2*abs(mu_x)./(mu_x^2+(fk-omega_x).^2);
exact_y = 2*abs(mu_y)./(mu_y^2+(fk-omega_y).^2);

%exact solutions for cross spectra
expression = 1./((LAMBDA(1))-sqrt(-1)*fk) + 1./(conj(LAMBDA(2))+sqrt(-1)*fk);
expression = expression*(-I1);


%% visualize

%plot power spectra
figure(1)
hold on
plot(fk,real(dad_x)/N*Delta,'-','color',[0.8500 0.3250 0.0980 .3],'linewidth',7)
plot(fk,real(dad_y)/N*Delta,'-','color',[0.4940 0.1840 0.5560 .3],'linewidth',7)
plot(fk,power_x/N*Delta,'-','color',[0.8500 0.3250 0.0980 .9],'linewidth',10)
plot(fk,exact_x,'-','color','y','linewidth',3)
plot(fk,power_y/N*Delta,'-','color',[0.4940 0.1840 0.5560 .9],'linewidth',10)
plot(fk,exact_y,'-','color','#FFC0CB','linewidth',3)
xlim([1 3.5])
ylim([0 20])
xlabel('frequency \nu')
ylabel('S_1(\nu)')
box on
axis square
set(gca,'fontsize',15)

%plot cross spectra
figure(2)
hold on
plot(fk,real(dad_cross)/N*Delta,'-','color',[.6 .6 .6 .6],'linewidth',10)
plot(fk,real(power_cross)/N*Delta,'-','color','#006400','linewidth',10)
plot(fk,real(expression),'-','color',[0 1 0],'linewidth',3)
xlim([1 3.5])
ylim([0 15])
xlabel('frequency \nu')
ylabel('Re($S_{\lambda_x,\lambda_y}$)', 'Interpreter', 'latex');
box on
axis square
set(gca,'fontsize',15)

figure(3)
hold on
plot(fk,imag(dad_cross)/N*Delta,'-','color',[.6 .6 .6 .6],'linewidth',10)
plot(fk,imag(power_cross)/N*Delta,'-','color',[0 0.4470 0.7410 .9],'linewidth',10)
plot(fk,imag(expression),'-','color',[0 1 1],'linewidth',3)
xlim([1 3.5])
ylim([-15 5])
xlabel('frequency \nu')
ylabel('Im($S_{\lambda_x,\lambda_y}$)', 'Interpreter', 'latex');
box on
axis square
set(gca,'fontsize',15)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% eigenvalue funcationsl

function[u] = compute_eig(gamma,D,w1,w2)

%cf stuff
M = 0;
N = 1;
cutoff = 250;

%initialize matrix
my_matrix_friend = zeros(2*cutoff+1,2*cutoff+1);

%form the matrix
count = 1;
for j=-cutoff+M:cutoff+M
    
    %diagonal entries
    my_matrix_friend(count,count) = Q_mid(N,j,D,w1,w2);
    
    %super/subdiagonal entries
    my_matrix_friend(count,count+1) = Q_plus(N,j,gamma);
    my_matrix_friend(count+1,count) = Q_minus(N,j+1,gamma);
    
    %update count
    count = count+1;
end

if mod(N,2) == 0
    my_matrix_friend = my_matrix_friend(1:2*cutoff+1,1:2*cutoff+1);
else
    my_matrix_friend(end,end) = Q_mid(N,cutoff+M+1,D,w1,w2);
end

%diagonalize
lambda_boi = eig(my_matrix_friend);

%sort
[~,temp] = sort(-real(lambda_boi));
lambda_boi = lambda_boi(temp);
lambda_boi = lambda_boi(1:2);
if abs(real(lambda_boi(1)) - real(lambda_boi(2))) < 1e-5 %before bifurcation
    if imag(lambda_boi(1)) > imag(lambda_boi(2))
        temp = lambda_boi(2);
        lambda_boi(2) = lambda_boi(1);
        lambda_boi(1) = temp;
    end
elseif abs(imag(lambda_boi(1)) - imag(lambda_boi(2))) < 1e-5 %after bifurcation
    if abs(real(lambda_boi(1))) > abs(real(lambda_boi(2)))
        temp = lambda_boi(2);
        lambda_boi(2) = lambda_boi(1);
        lambda_boi(1) = temp;
    end
end
u = lambda_boi(1:2);
end


%% eigenfunction funcational

function[cf,cb] = compute_efunc(M,gamma,D,w1,w2,lambda)

%cf stufuf
N = 1;
cutoff = 50;
S = zeros(1,cutoff-1);
ccount = 1;

%loop
for n=M:cutoff+M
    
    %initialize
    a = zeros(1,cutoff);
    b = zeros(1,cutoff);
    
    %initial values follow a different pattern than the rest
    a(1) = -Q_minus(N,n+1,gamma);
    
    %create b vector
    count = 1;
    for j=n:cutoff+n
        b(count) = -lambda + Q_mid(N,j+1,D,w1,w2);
        count = count+1;
    end
    
    %create a vector
    count = 2;
    for j=n:cutoff+n
        a(count) = -Q_plus(N,j+1,gamma)*Q_minus(N,j+2,gamma);
        count = count+1;
    end
    
    %initialize solution vectors
    A = zeros(1,cutoff);
    B = zeros(1,cutoff);
    
    %initial conditions
    A(1) = 1;
    A(2) = 0;
    B(1) = 0;
    B(2) = 1;
    
    %loop
    for j=1:cutoff-2
        A(j+2) = b(j)*A(j+1) + a(j)*A(j);
        B(j+2) = b(j)*B(j+1) + a(j)*B(j);
    end
    
    %output Nth approximant
    S(ccount) = A(end)/B(end);
    
    %indexing
    ccount = ccount + 1;
end

%coefficient vector
cf = zeros(1,cutoff);
cf(1) = 1;

%solve for coefficients
for j=1:cutoff
    cf(j+1) = S(j)*cf(j);
end

%reset
S = 0;

%indexing L :(
ccount = 1;

%loop
for n=M:-1:-cutoff+M
    
    %initialize
    a = zeros(1,cutoff);
    b = zeros(1,cutoff);
    
    %initial values follow a different pattern than the rest
    a(1) = -Q_plus(N,n-1,gamma);
    
    %create b vector
    count = 1;
    for j=n:-1:-cutoff+n
        b(count) = -lambda + Q_mid(N,j-1,D,w1,w2);
        count = count+1;
    end
    
    %create a vector
    count = 2;
    for j=n:-1:-cutoff+n
        a(count) = -Q_minus(N,j-1,gamma)*Q_plus(N,j-2,gamma);
        count = count+1;
    end
    
    %initialize solution vectors
    A = zeros(1,cutoff);
    B = zeros(1,cutoff);
    
    %initial conditions
    A(1) = 1;
    A(2) = 0;
    B(1) = 0;
    B(2) = 1;
    
    %loop
    for j=1:cutoff-2
        A(j+2) = b(j)*A(j+1) + a(j)*A(j);
        B(j+2) = b(j)*B(j+1) + a(j)*B(j);
    end
    
    %output Nth approximant
    S(ccount) = A(end)/B(end);
    
    %update indexing
    ccount = ccount + 1;
    
end

%coefficient vector
cb = zeros(1,cutoff);
cb(1) = 1;

%solve for coefficients
for j=1:cutoff
    cb(j+1) = S(j)*cb(j);
end

end


%% compute Q function at specified points

function[Q] = computeQ(x,y,cf,cb,M)

%construct
cutoff = 50;
N = 1;
Qf = 0;
count = 1;
for n=M:cutoff+M
    Qf = Qf + cf(count).*exp(sqrt(-1)*(n*x+(N-n)*y));
    count = count+1;
end
Qb = 0;
count = 2;
for n=M-1:-1:-cutoff+M
    Qb = Qb + cb(count).*exp(sqrt(-1)*(n*x+(N-n)*y));
    count = count + 1;
end
Q = (Qf+Qb);


end


%% eigenfunction funcational

function[cf,cb] = compute_efunc_forward(M,gamma,D,w1,w2,lambda)

%cf stufuf
N = 0;
cutoff = 50;
S = zeros(1,cutoff-1);
ccount = 1;

%loop
for n=M:cutoff+M
    
    %initialize
    a = zeros(1,cutoff);
    b = zeros(1,cutoff);
    
    %initial values follow a different pattern than the rest
    a(1) = -Q_minus_forward(N,n+1,gamma);
    
    %create b vector
    count = 1;
    for j=n:cutoff+n
        b(count) = -lambda + Q_mid_forward(N,j+1,D,w1,w2);
        count = count+1;
    end
    
    %create a vector
    count = 2;
    for j=n:cutoff+n
        a(count) = -Q_plus_forward(N,j+1,gamma)*Q_minus_forward(N,j+2,gamma);
        count = count+1;
    end
    
    %initialize solution vectors
    A = zeros(1,cutoff);
    B = zeros(1,cutoff);
    
    %initial conditions
    A(1) = 1;
    A(2) = 0;
    B(1) = 0;
    B(2) = 1;
    
    %loop
    for j=1:cutoff-2
        A(j+2) = b(j)*A(j+1) + a(j)*A(j);
        B(j+2) = b(j)*B(j+1) + a(j)*B(j);
    end
    
    %output Nth approximant
    S(ccount) = A(end)/B(end);
    
    %indexing
    ccount = ccount + 1;
end

%coefficient vector
cf = zeros(1,cutoff);
cf(1) = 1;

%solve for coefficients
for j=1:cutoff
    cf(j+1) = S(j)*cf(j);
end

%reset
S = 0;

%indexing L :(
ccount = 1;

%loop
for n=M:-1:-cutoff+M
    
    %initialize
    a = zeros(1,cutoff);
    b = zeros(1,cutoff);
    
    %initial values follow a different pattern than the rest
    a(1) = -Q_plus_forward(N,n-1,gamma);
    
    %create b vector
    count = 1;
    for j=n:-1:-cutoff+n
        b(count) = -lambda + Q_mid_forward(N,j-1,D,w1,w2);
        count = count+1;
    end
    
    %create a vector
    count = 2;
    for j=n:-1:-cutoff+n
        a(count) = -Q_minus_forward(N,j-1,gamma)*Q_plus_forward(N,j-2,gamma);
        count = count+1;
    end
    
    %initialize solution vectors
    A = zeros(1,cutoff);
    B = zeros(1,cutoff);
    
    %initial conditions
    A(1) = 1;
    A(2) = 0;
    B(1) = 0;
    B(2) = 1;
    
    %loop
    for j=1:cutoff-2
        A(j+2) = b(j)*A(j+1) + a(j)*A(j);
        B(j+2) = b(j)*B(j+1) + a(j)*B(j);
    end
    
    %output Nth approximant
    S(ccount) = A(end)/B(end);
    
    %update indexing
    ccount = ccount + 1;
    
end

%coefficient vector
cb = zeros(1,cutoff);
cb(1) = 1;

%solve for coefficients
for j=1:cutoff
    cb(j+1) = S(j)*cb(j);
end

end


%% compute P0 function at specified points

function[Q] = computeP0(x,y,cf,cb,M)

%construct
cutoff = 50;
N = 0;
Qf = 0;
count = 1;
for n=M:cutoff+M
    Qf = Qf + cf(count).*exp(sqrt(-1)*(n*x+(N-n)*y));
    count = count+1;
end
Qb = 0;
count = 2;
for n=M-1:-1:-cutoff+M
    Qb = Qb + cb(count).*exp(sqrt(-1)*(n*x+(N-n)*y));
    count = count + 1;
end
Q = (Qf+Qb);


end


%% cf functional

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


%% cf (forward) functional

%Q-
function[u] = Q_minus_forward(N,m,gamma)
u = gamma/2*(2*m-N);
end

%Q
function[u] = Q_mid_forward(N,m,D,w1,w2)
u = -D*(m^2+(N-m)^2)-sqrt(-1)*w1*m - sqrt(-1)*w2*(N-m);
end

%Q+
function[u] = Q_plus_forward(N,m,gamma)
u = gamma/2*(N-2*m);
end


