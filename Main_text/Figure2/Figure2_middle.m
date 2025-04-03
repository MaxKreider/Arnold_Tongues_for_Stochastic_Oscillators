%% produce isochrons for both and eigenvalues for 1 of the oscillators

%keep it clean
clc
clf

%specify rectangular domain
a = 0;
b = 2*pi;
c = 0;
d = 2*pi;

%specify (full) grid size
N = 400+1;
M = 400+1;

%mesh
x = linspace(a,b,N);
y = linspace(c,d,M);
[X,Y] = meshgrid(x,y);

%Q function
Q = cos(X) - sqrt(-1)*sin(X) + 0*Y;
Q = Q*exp(sqrt(-1)*pi);

%isochrons and isostable function
fig = figure(1);
Qtest = reshape(Q,[N, M])';
contour(X,Y,angle(Qtest)+pi,13,'LineWidth',3)
colorbar
colormap jet
xlabel('time t')
ylabel('x(t), y(t)')
box on
set(gca,'FontSize',15)
xlim([0 2*pi])
ylim([0 2*pi])
axis square


%% generate phase plane and trajectories

for pq=1:30
    
    %parameters
    D = .1;
    w1 = 2;
    w2 = 2.5;
    
    %time
    Delta = 1/500;
    N = 2^11+1000;
    t = 0:Delta:(N-1)*Delta;
    
    %solution vectors
    XXt = zeros(1,length(t));
    YYt = zeros(1,length(t));
    
    %trajectory initial value (uniformly random on the torus)
    SOL = [0;0];
    
    %initial perturbed phase
    XXt(1) = SOL(1);
    YYt(1) = SOL(2);
    
    %compute phases
    for pp=1:length(t)-1
        
        %EM method
        SOL = SOL + Delta*[w1;w2] + sqrt(Delta)*sqrt(2*D)*randn(2,1);
        
        %compute Q function at that point
        XXt(pp+1) = SOL(1);
        YYt(pp+1) = SOL(2);
        
    end
    
    %phase plane
    figure(1);
    hold on
    scatter(mod(t,2*pi), mod(XXt,2*pi), '.', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
    scatter(mod(t,2*pi), mod(YYt,2*pi), '.', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7], 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
    
end

%time
Delta = 1/500;
N = 2^16;
t = 0:Delta:(N-1)*Delta;

%solution vectors
XXt = zeros(1,length(t));
YYt = zeros(1,length(t));

%trajectory initial value (uniformly random on the torus)
SOL = [0;0];

%initial perturbed phase
XXt(1) = SOL(1);
YYt(1) = SOL(2);

%compute phases
for pp=1:length(t)-1
    
    %EM method
    SOL = SOL + Delta*[w1;w2] + sqrt(Delta)*sqrt(2*D)*randn(2,1);
    
    %compute Q function at that point
    XXt(pp+1) = SOL(1);
    YYt(pp+1) = SOL(2);
    
end
    
%time series
figure(2);
subplot(2,1,2)
hold on
scatter(t, mod(XXt,2*pi), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
scatter(t, mod(YYt,2*pi), 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7], 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
xlabel('time t')
ylabel('x(t)')
set(gca,'fontsize',12)
xlim([0 33])
set(gca,'FontSize',15)


%% eigenvalues

%uncoupled eigenvalues
lambda_x = zeros(5,5);
for i=-2:2
    for j=-2:2
        lambda_x(i+3,j+3) = -.1*(i^2+0*j^2) + sqrt(-1)*w1*(i+0*j);
    end
end
lambda_y = zeros(5,5);
for i=-2:2
    for j=-2:2
        lambda_y(i+3,j+3) = -.1*(0*i^2+j^2) + sqrt(-1)*w2*(0*i+j);
    end
end

%eigenvalue plot
figure(2)
subplot(2,1,1)
hold on
grid on
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
box on
set(gca,'FontSize',15)
xlim([-1 0])
plot(real(lambda_x),imag(lambda_x),'.','color','k','MarkerSize',40)
plot(real(lambda_y),imag(lambda_y),'.','color',[.7 .7 .7],'MarkerSize',35)


%% compute numerical trajectories

%system parameters
w1 = 2;
w2 = 2.5;
D = 0.1;

%variable system parameters
gamma = 0;

%setup for psd do-da
Delta = 1/100;
N = 2^16;   
t = 0:Delta:(N-1)*Delta;

%frequency vector
k = (-N/2:N/2-1);
fk = 1/(N*Delta)*k*2*pi;

%number of trials
M = 2000;


%% do it

%initialize
dad_x = 0;
dad_y = 0;
    
%% compute numerical trajectories

parfor qq=1:M
    
    %solution vector
    XX = zeros(1,length(t));
    YY = zeros(1,length(t));
    
    %trajectory initial value (uniformly random on the torus)
    SOL = rand(2,1)*2*pi;
    
    %initial perturbed phase
    XX(1) = cos(SOL(1)) + sqrt(-1)*sin(SOL(1));
    YY(1) = cos(SOL(2)) + sqrt(-1)*sin(SOL(2));
    
    %compute phases
    for pp=1:length(t)-1
        
        %EM method
        SOL = SOL + Delta*[w1+gamma*sin(SOL(2)-SOL(1));w2+gamma*sin(SOL(1)-SOL(2))] + sqrt(Delta)*sqrt(2*D)*randn(2,1);
        
        %compute Q function at that point
        XX(pp+1) = cos(SOL(1)) + sqrt(-1)*sin(SOL(1));
        YY(pp+1) = cos(SOL(2)) + sqrt(-1)*sin(SOL(2));
        
    end
    
    %take fft
    Fc_X = fftshift(fft(XX));
    Fc_Y = fftshift(fft(YY));
    
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
semilogy(fk,real(dad_x)/N*Delta,'-','color','k','linewidth',2)
semilogy(fk,real(dad_y)/N*Delta,'-','color',[.7 .7 .7],'linewidth',2)
xlim([1 3.5])
ylim([10^-2 10^2])
xlabel('frequency \nu')
ylabel('S_1(\nu)')
box on
axis square
set(gca,'FontSize',15)
set(gca, 'YScale', 'log')
