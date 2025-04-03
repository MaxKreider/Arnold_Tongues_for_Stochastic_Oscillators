%% produce isochrons for both and eigenvalues for 1 of the oscillators

%keep it clean
clc
clf

%specify rectangular domain
a = -4;
b = 4;
c = -4;
d = 4;

%specify (full) grid size
N = 400+1;
M = 400+1;

%mesh
x = linspace(a,b,N);
y = linspace(c,d,M);
[X,Y] = meshgrid(x,y);

%step size
h = (b-a)/(N-1);
k = (d-c)/(M-1);

%parameters
D = .1;
K = .1;
w = 2;

%variable coefficients
f = @(x,y) D + 0*x.*y;
g = @(x,y) D + 0*x.*y;
m = @(x,y) -K*x + w*y + 0*x.*y;
n = @(x,y) -w*x + -K*y + 0*x.*y;


%% discretize operator

%set up variable coefficients
X1 = reshape(X(1:end,1:end)',N*M,1);
Y1 = reshape(Y(1:end,1:end)',N*M,1);

%specify variable coefficient values
V1 = f(X1,Y1);
V2 = g(X1,Y1);
V3 = m(X1,Y1);
V4 = n(X1,Y1);

%specify coefficient vectors
coeff_A = 2*(V1+V2*(h/k)^2);
coeff_B = -V1 + 1/2*h*V3;
coeff_C = -V1 - 1/2*h*V3;
coeff_D = -(h/k)^2*V2 - h^2/(2*k)*V4;
coeff_E = -(h/k)^2*V2 + h^2/(2*k)*V4;

coeff_a = 1/12*V1 - h/12*V3;
coeff_b = -16/12*V1 + h*8/12*V3;
coeff_c = -16/12*V1 - h*8/12*V3;
coeff_d = 1/12*V1 + h/12*V3;
coeff_i = 1/12*(h/k)^2*V2 - 1/12*(h^2/k)*V4;
coeff_h = -16/12*(h/k)^2*V2 + 8/12*(h^2/k)*V4;
coeff_f = -16/12*(h/k)^2*V2 - 8/12*(h^2/k)*V4;
coeff_e = 1/12*(h/k)^2*V2 + 1/12*(h^2/k)*V4;
coeff_j = 30/12*(V1+(h/k)^2*V2);

%modify coeff_D for neumann cond
coeff_D(1:N) = coeff_D(1:N) + coeff_E(1:N);

%modify coeff_E for neumann cond
coeff_E(end-(N-1):end) = coeff_E(end-(N-1):end) + coeff_D(end-(N-1):end);

%modify coeff_B for neumann cond
coeff_B(N:N:end) = coeff_B(N:N:end) + coeff_C(N:N:end);

%modify coeff_C for neumann and hope for singular
coeff_C(1:N:end) = coeff_C(1:N:end) + coeff_B(1:N:end);

%construct outer A
e = ones((N)*(M),1);
e(2*(N)+1:(M-2)*(N),1:end)=0;
e(2*(N)+1:(N):(M-2)*(N),1:end)=1;
e(3*(N):(N):(M-2)*(N),1:end)=1;
e(2*(N)+2:(N):(M-2)*(N),1:end)=1;
e(3*(N)-1:(N):(M-2)*(N),1:end)=1;
e2 = e;
e3 = e;
e2(N:N:end) = 0;
e3(N+1:N:end) = 0;

A_out = spdiags([coeff_D.*e coeff_C.*e2 coeff_A.*e coeff_B.*e3 coeff_E.*e],[-(N) -1 0 1 (N)],(N)*(M),(N)*(M))'; 

%construct inner A
e = zeros((N)*(M),1);
e(2*(N)+3:(M-2)*(N)-2,1:end)=1;
e(2*(N)+1:(N):(M-2)*(N),1:end)=0;
e(3*(N):(N):(M-2)*(N),1:end)=0;
e(2*(N)+2:(N):(M-2)*(N),1:end)=0;
e(3*(N)-1:(N):(M-2)*(N),1:end)=0;
e2 = e;
e3 = e;
e4 = e;
e5 = e;
e6 = e;
e7 = e;
e2(N:N:end) = 0;
e3(N+1:N:end) = 0;
e4(N-1:N:end) = 0;
e5(N+2:N:end) = 0;

A_in = spdiags([coeff_e.*e coeff_f.*e coeff_d.*e4 coeff_c.*e2 coeff_j.*e coeff_b.*e3 coeff_a.*e5 coeff_h.*e coeff_i.*e],[-2*(N) -(N) -2 -1 0 1 2 (N) 2*(N)],(N)*(M),(N)*(M))'; 

%construct backwards A
A_b = -1/h^2*(A_in + A_out);

%forward matrix
A_f = A_b.';
            
%11 smallest evalues and time it
tic
[VV_b, lambda_b] = eigs(A_b,12,1e-4);
[VV_f, lambda_f] = eigs(A_f,12,1e-4);
toc

%eigenvalue plot
figure(2)
subplot(2,1,1)
lambda_b=diag(lambda_b)
hold on
plot(real(lambda_b),imag(lambda_b),'.','color','#006400','MarkerSize',40)
plot(real(lambda_f),imag(lambda_f),'.','color','#006400','MarkerSize',40)
grid on
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
box on
set(gca,'FontSize',15)
xlim([-1 0])

%unnormalize to plot
Q_b = VV_b./max(VV_b, [], 1);
Q_f = VV_f./max(VV_f, [], 1);

%isochrons and isostable function
fig = figure(1);
Qtest = reshape(Q_b(:,7),[N, M])';
contour(X,Y,angle(Qtest)+pi,13,'LineWidth',3)
colorbar
colormap jet
xlabel('x')
ylabel('y')
box on
set(gca,'FontSize',12)
xlim([-3 3])
ylim([-3 3])
axis square
set(gca,'FontSize',15)


%% generate phase plane and trajectories

%parameters
D = .1;
k = .1;
w1 = 2;
w2 = 2.5;
gamma = 0;

%construct A
A = [-k-gamma, w1, gamma, 0;
    -w1, -k-gamma, 0, gamma;
    gamma, 0, -k-gamma, w2;
    0, gamma, -w2, -k-gamma];

%time
Delta = 1/500;
N = 2^15;
t = 0:Delta:(N-1)*Delta;

%solution vectors
XXt = zeros(1,length(t));
YYt = zeros(1,length(t));
XXt2 = zeros(1,length(t));
YYt2 = zeros(1,length(t));

%trajectory initial value (uniformly random on the torus)
SOL = [0;0;0;0];

%initial perturbed phase
XXt(1) = SOL(1);
YYt(1) = SOL(3);
XXt2(1) = SOL(2);
YYt2(1) = SOL(4);

%compute phases
for pp=1:length(t)-1
    
    %EM method
    SOL = SOL + Delta*A*SOL + sqrt(Delta)*sqrt(2*D)*randn(4,1);
    
    %compute Q function at that point
    XXt(pp+1) = SOL(1);
    YYt(pp+1) = SOL(3);
    XXt2(pp+1) = SOL(2);
    YYt2(pp+1) = SOL(4);
    
end

%phase plane
figure(1);
hold on
plot(XXt,XXt2,'-','color','#006400')
plot(YYt,YYt2,'-','color','g')

%time series
figure(2);
subplot(2,1,2)
hold on
plot(t(1:round(end/2)),XXt(1:round(end/2)),'-','color','#006400','linewidth',2)
plot(t(1:round(end/2)),YYt(1:round(end/2)),'-','color','g','linewidth',2)
xlabel('time t')
ylabel('x(t)')
set(gca,'fontsize',15)
xlim([0 33])


%% generate other eigenvalues

%specify rectangular domain
a = -4;
b = 4;
c = -4;
d = 4;

%specify (full) grid size
N = 400+1;
M = 400+1;

%mesh
x = linspace(a,b,N);
y = linspace(c,d,M);
[X,Y] = meshgrid(x,y);

%step size
h = (b-a)/(N-1);
k = (d-c)/(M-1);

%parameters
D = .1;
K = .1;
w = 2.5;

%variable coefficients
f = @(x,y) D + 0*x.*y;
g = @(x,y) D + 0*x.*y;
m = @(x,y) -K*x + w*y + 0*x.*y;
n = @(x,y) -w*x + -K*y + 0*x.*y;


%% discretize operator

%set up variable coefficients
X1 = reshape(X(1:end,1:end)',N*M,1);
Y1 = reshape(Y(1:end,1:end)',N*M,1);

%specify variable coefficient values
V1 = f(X1,Y1);
V2 = g(X1,Y1);
V3 = m(X1,Y1);
V4 = n(X1,Y1);

%specify coefficient vectors
coeff_A = 2*(V1+V2*(h/k)^2);
coeff_B = -V1 + 1/2*h*V3;
coeff_C = -V1 - 1/2*h*V3;
coeff_D = -(h/k)^2*V2 - h^2/(2*k)*V4;
coeff_E = -(h/k)^2*V2 + h^2/(2*k)*V4;

coeff_a = 1/12*V1 - h/12*V3;
coeff_b = -16/12*V1 + h*8/12*V3;
coeff_c = -16/12*V1 - h*8/12*V3;
coeff_d = 1/12*V1 + h/12*V3;
coeff_i = 1/12*(h/k)^2*V2 - 1/12*(h^2/k)*V4;
coeff_h = -16/12*(h/k)^2*V2 + 8/12*(h^2/k)*V4;
coeff_f = -16/12*(h/k)^2*V2 - 8/12*(h^2/k)*V4;
coeff_e = 1/12*(h/k)^2*V2 + 1/12*(h^2/k)*V4;
coeff_j = 30/12*(V1+(h/k)^2*V2);

%modify coeff_D for neumann cond
coeff_D(1:N) = coeff_D(1:N) + coeff_E(1:N);

%modify coeff_E for neumann cond
coeff_E(end-(N-1):end) = coeff_E(end-(N-1):end) + coeff_D(end-(N-1):end);

%modify coeff_B for neumann cond
coeff_B(N:N:end) = coeff_B(N:N:end) + coeff_C(N:N:end);

%modify coeff_C for neumann and hope for singular
coeff_C(1:N:end) = coeff_C(1:N:end) + coeff_B(1:N:end);

%construct outer A
e = ones((N)*(M),1);
e(2*(N)+1:(M-2)*(N),1:end)=0;
e(2*(N)+1:(N):(M-2)*(N),1:end)=1;
e(3*(N):(N):(M-2)*(N),1:end)=1;
e(2*(N)+2:(N):(M-2)*(N),1:end)=1;
e(3*(N)-1:(N):(M-2)*(N),1:end)=1;
e2 = e;
e3 = e;
e2(N:N:end) = 0;
e3(N+1:N:end) = 0;

A_out = spdiags([coeff_D.*e coeff_C.*e2 coeff_A.*e coeff_B.*e3 coeff_E.*e],[-(N) -1 0 1 (N)],(N)*(M),(N)*(M))'; 

%construct inner A
e = zeros((N)*(M),1);
e(2*(N)+3:(M-2)*(N)-2,1:end)=1;
e(2*(N)+1:(N):(M-2)*(N),1:end)=0;
e(3*(N):(N):(M-2)*(N),1:end)=0;
e(2*(N)+2:(N):(M-2)*(N),1:end)=0;
e(3*(N)-1:(N):(M-2)*(N),1:end)=0;
e2 = e;
e3 = e;
e4 = e;
e5 = e;
e6 = e;
e7 = e;
e2(N:N:end) = 0;
e3(N+1:N:end) = 0;
e4(N-1:N:end) = 0;
e5(N+2:N:end) = 0;

A_in = spdiags([coeff_e.*e coeff_f.*e coeff_d.*e4 coeff_c.*e2 coeff_j.*e coeff_b.*e3 coeff_a.*e5 coeff_h.*e coeff_i.*e],[-2*(N) -(N) -2 -1 0 1 2 (N) 2*(N)],(N)*(M),(N)*(M))'; 

%construct backwards A
A_b = -1/h^2*(A_in + A_out);

%forward matrix
A_f = A_b.';
            
%11 smallest evalues and time it
tic
[VV_b, lambda_b] = eigs(A_b,13,1e-4);
[VV_f, lambda_f] = eigs(A_f,13,1e-4);
toc

%eigenvalue plot
figure(2)
subplot(2,1,1)
lambda_b=diag(lambda_b)
hold on
plot(real(lambda_b),imag(lambda_b),'.','color','g','MarkerSize',35)
plot(real(lambda_f),imag(lambda_f),'.','color','g','MarkerSize',35)


%% compute numerical trajectories

%stuff
M = 5000;
dad_x = 0;
dad_y = 0;

%setup for psd do-da
Delta = 1/100;
N = 2^18;
t = 0:Delta:(N-1)*Delta;

%frequency vector
kkk = (-N/2:N/2-1);
fk = 1/(N*Delta)*kkk*2*pi;

parfor qq=1:M
    
    %solution vectors
    XXt = zeros(1,length(t));
    YYt = zeros(1,length(t));
    
    %trajectory initial value (uniformly random on the torus)
    SOL = randn(4,1);
    
    %initial perturbed phase
    XXt(1) = SOL(1);
    YYt(1) = SOL(3);
    
    %compute phases
    for pp=1:length(t)-1
        
        %EM method
        SOL = SOL + Delta*A*SOL + sqrt(Delta)*sqrt(2*D)*randn(4,1);
        
        %compute Q function at that point
        XXt(pp+1) = SOL(1);
        YYt(pp+1) = SOL(3);
        
    end
    
    %take fft
    Fc_X = fftshift(fft(XXt));
    Fc_Y = fftshift(fft(YYt));
    
    %power spectra
    dad_x = dad_x + abs(Fc_X).^2;
    dad_y = dad_y + abs(Fc_Y).^2;
    
end

%normalize
dad_x = dad_x/M;
dad_y = dad_y/M;

%plot power spectra
figure(3)
hold on
semilogy(fk,real(dad_x)/N*Delta,'-','color','#006400','linewidth',2)
semilogy(fk,real(dad_y)/N*Delta,'-','color','g','linewidth',2)
xlim([1 3.5])
ylim([10^-2 10^2])
xlabel('frequency \nu')
ylabel('S_1(\nu)')
box on
axis square
set(gca,'fontsize',15)
set(gca, 'YScale', 'log')

