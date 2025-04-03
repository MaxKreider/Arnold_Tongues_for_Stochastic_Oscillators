%% first osc evals

%mr clean
clc
clf
clear

%size
SIZE = 3;

%stoichiometric matrix
S = -eye(SIZE)*1 + 1*circshift(eye(SIZE),-1);

%initial amount of each species
x0 = zeros(SIZE,1);
x0(1) = 1;

%define reaction rates
a = 1;

%interpolation time vector and Tmax
Delta = 1/5000;
N = 2^15;
t_query = 0:Delta:(N-1)*Delta;
Tmax = t_query(end) + 1000*Delta;

%number of trials
bigN = 1;

%frequency vector
kkk = (-N/2:N/2-1);
fk = 1/(N*Delta)*kkk*2*pi;

%number of repititions for the power spectra
M = 1000;

%yes
power_x = 0;
power_x_Q = 0;


%form the matrix
A = (-eye(SIZE)*a + a*circshift(eye(SIZE),-1))';
[v,lambda] = eig(A);
[w,lambda] = eig(A.');
lambda = diag(lambda)
imag(lambda)./real(lambda)
W = w(:,1)
P0 = v(:,3)/sum(v(:,3))

%plot the spectrum
figure(69)
plot(real(lambda),imag(lambda),'m.','markersize',40)
hold on
grid on
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
box on
set(gca,'FontSize',15)
xlim([-2 0])


%% second osc evals

%size
SIZE = 3;

%stoichiometric matrix
S = -eye(SIZE)*1 + 1*circshift(eye(SIZE),-1);

%initial amount of each species
x0 = zeros(SIZE,1);
x0(1) = 1;

%define reaction rates
a = 1.5;

%interpolation time vector and Tmax
Delta = 1/5000;
N = 2^15;
t_query = 0:Delta:(N-1)*Delta;
Tmax = t_query(end) + 1000*Delta;

%number of trials
bigN = 1;

%frequency vector
kkk = (-N/2:N/2-1);
fk = 1/(N*Delta)*kkk*2*pi;

%number of repititions for the power spectra
M = 1000;

%yes
power_x = 0;
power_x_Q = 0;


%form the matrix
A = (-eye(SIZE)*a + a*circshift(eye(SIZE),-1))';
[v,lambda] = eig(A);
[w,lambda] = eig(A.');
lambda = diag(lambda)
imag(lambda)./real(lambda)
W = w(:,1)
P0 = v(:,3)/sum(v(:,3))

%plot the spectrum
figure(69)
plot(real(lambda),imag(lambda),'.','color',[255, 182, 193 256]/256,'markersize',35)
hold on
grid on
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
box on
set(gca,'FontSize',15)
xlim([-4 0])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% first osc sims

%size
SIZE = 3;

%stoichiometric matrix
S = -eye(SIZE)*1 + 1*circshift(eye(SIZE),-1);

%initial amount of each species
x0 = zeros(SIZE,1);
x0(1) = 1;

%define reaction rates
a = 1;

%interpolation time vector and Tmax
Delta = 1/500;
N = 2^15;
t_query = 0:Delta:(N-1)*Delta;
Tmax = t_query(end) + 1000*Delta;

%number of trials
bigN = 1;

%frequency vector
kkk = (-N/2:N/2-1);
fk = 1/(N*Delta)*kkk*2*pi;

%yes
power_x = 0;
power_x_Q = 0;

%form the matrix
A = (-eye(SIZE)*a + a*circshift(eye(SIZE),-1))';
[v,lambda] = eig(A);
[w,lambda] = eig(A.');
lambda = diag(lambda)
imag(lambda)./real(lambda)
W = w(:,1)
P0 = v(:,3)/sum(v(:,3))

%solution vector
SOLBOI = zeros(SIZE,length(t_query));
SOLQ = zeros(1,length(t_query));

%loopin it
for i=1:bigN
    
    %define propensity function
    prop = @(x) a*[x(1);x(2);x(3)];
    
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
        alpha=prop(X(:,count));
        xi=rand;
        s=-1/sum(alpha)*log(1-xi);
        
        %draw the reaction
        phi=0;
        ell=0;
        eta=rand;
        while (phi < eta)
            ell=ell+1;
            phi=phi + alpha(ell)/sum(alpha);
        end
        
        %update
        t=t+s;
        X(:,count+1)=X(:,count) + S(:,ell);
        tspan(count+1)=t;
        
        %update count
        count = count + 1;
    end
    
    %do the smush smush
    X_sol = zeros(1,length(X));
    for pp=1:length(X_sol)
        X_sol(pp) = find(X(:,pp)==1);
    end
    
end

%visualize
figure(1)
hold on
stairs(tspan,flip(X_sol),'m','LineWidth',3)
ylim([-1 4])
xlim([0 15])


%% first osc sims

%size
SIZE = 3;

%stoichiometric matrix
S = -eye(SIZE)*1 + 1*circshift(eye(SIZE),-1);

%initial amount of each species
x0 = zeros(SIZE,1);
x0(1) = 1;

%define reaction rates
a = 1.5;

%interpolation time vector and Tmax
Delta = 1/500;
N = 2^15;
t_query = 0:Delta:(N-1)*Delta;
Tmax = t_query(end) + 1000*Delta;

%number of trials
bigN = 1;

%frequency vector
kkk = (-N/2:N/2-1);
fk = 1/(N*Delta)*kkk*2*pi;

%yes
power_x = 0;
power_x_Q = 0;

%form the matrix
A = (-eye(SIZE)*a + a*circshift(eye(SIZE),-1))';
[v,lambda] = eig(A);
[w,lambda] = eig(A.');
lambda = diag(lambda)
imag(lambda)./real(lambda)
W = w(:,1)
P0 = v(:,3)/sum(v(:,3))

%solution vector
SOLBOI = zeros(SIZE,length(t_query));
SOLQ = zeros(1,length(t_query));

%loopin it
for i=1:bigN
    
    %define propensity function
    prop = @(x) a*[x(1);x(2);x(3)];
    
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
        alpha=prop(X(:,count));
        xi=rand;
        s=-1/sum(alpha)*log(1-xi);
        
        %draw the reaction
        phi=0;
        ell=0;
        eta=rand;
        while (phi < eta)
            ell=ell+1;
            phi=phi + alpha(ell)/sum(alpha);
        end
        
        %update
        t=t+s;
        X(:,count+1)=X(:,count) + S(:,ell);
        tspan(count+1)=t;
        
        %update count
        count = count + 1;
    end
    
    %do the smush smush
    X_sol = zeros(1,length(X));
    for pp=1:length(X_sol)
        X_sol(pp) = find(X(:,pp)==1);
    end
    
end

%visualize
figure(1)
hold on
stairs(tspan,flip(X_sol),'color',[255, 182, 193 256]/256,'LineWidth',2.5)
ylim([0 4])
xlim([0 15])
xlabel('time t')
ylabel('X_A, X_B')
set(gca,'FontSize',15)
box on
axis square


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% first osc power spectra

%size
SIZE = 3;

%stoichiometric matrix
S = -eye(SIZE)*1 + 1*circshift(eye(SIZE),-1);

%initial amount of each species
x0 = zeros(SIZE,1);
x0(1) = 1;

%define reaction rates
a = 1;

%interpolation time vector and Tmax
Delta = 1/5000;
N = 2^15;
t_query = 0:Delta:(N-1)*Delta;
Tmax = t_query(end) + 10000*Delta;

%number of trials
bigN = 1;

%frequency vector
kkk = (-N/2:N/2-1);
fk = 1/(N*Delta)*kkk*2*pi;

%number of repititions for the power spectra
M = 200;

%yes
power_x = 0;
power_x_Q = 0;

%form the matrix
A = (-eye(SIZE)*a + a*circshift(eye(SIZE),-1))';
[v,lambda] = eig(A);
[w,lambda] = eig(A.');
lambda = diag(lambda)
imag(lambda)./real(lambda)
W = w(:,1)
P0 = v(:,3)/sum(v(:,3))

lambda = lambda(1);


%loop for power spectra
for mmm = 1:M
    
    %solution vector
    SOLBOI = zeros(SIZE,length(t_query));
    SOLQ = zeros(1,length(t_query));
    
    %loopin it
    for i=1:bigN
        
        %define propensity function
        prop = @(x) a*[x(1);x(2);x(3)];
        
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
            alpha=prop(X(:,count));
            xi=rand;
            s=-1/sum(alpha)*log(1-xi);
            
            %draw the reaction
            phi=0;
            ell=0;
            eta=rand;
            while (phi < eta)
                ell=ell+1;
                phi=phi + alpha(ell)/sum(alpha);
            end
            
            %update
            t=t+s;
            X(:,count+1)=X(:,count) + S(:,ell);
            tspan(count+1)=t;
            
            %update count
            count = count + 1;
        end
        
        
        %do the smush smush
        X_sol = zeros(1,length(X));
        for pp=1:length(X_sol)
            X_sol(pp) = find(X(:,pp)==1);
        end
        
        %interpolate the data onto an evenly spaced time grid
        
        X_interped = interp1(tspan,X_sol,t_query,'previous');
        
       
        
        %do it
        for pp = 1:length(t_query)
            idx = X_interped(pp);
            SOLBOI(idx, pp) = SOLBOI(idx, pp) + 1;
        end
        
    end
    
    %normalize the solution to reflect the percentage of realizations in a
    %state
    SOLBOI = SOLBOI/bigN;
    

    %% power spectra
    
    %do the Q thing
    for yy=1:length(t_query)
        SOLQ(yy) = dot(W,SOLBOI(:,yy));
        SOLQ(yy) = SOLQ(yy)/(sum(abs(SOLQ(yy))));
    end
    
    %fft
    Fc_x = fftshift(fft(SOLBOI(2,:) - 0*mean(SOLBOI(2,:))));
    Fc_Qx = fftshift(fft(SOLQ - 0*mean(SOLQ)));
    
    %power spectra
    power_x = power_x + abs(Fc_x).^2;
    power_x_Q = power_x_Q + abs(Fc_Qx).^2;
end

%normalize those rascals
power_x = power_x/M;
power_x_Q = power_x_Q/M;

%exact
exact = -2*real(lambda)./(real(lambda)^2 + (fk-imag(lambda)).^2);
maxmax = max(exact);

%normalize
power_x_Q = power_x_Q/N*Delta;
%pdaddy_x_Q = pdaddy_x_Q/max(pdaddy_x_Q)*maxmax;

power_x = power_x/N*Delta;
%pdaddy_x = pdaddy_x/max(pdaddy_x)*maxmax;
power_x(power_x>0.4) = nan;

%plot power spectra
figure(7346)
hold on
semilogy(fk,(power_x),'-','color','m','linewidth',2)
xlabel('frequency \nu')
ylabel('S_1(\nu)')
box on
axis square
set(gca,'FontSize',15)
xlim([-150 150])


%% second osc power spectra

%size
SIZE = 3;

%stoichiometric matrix
S = -eye(SIZE)*1 + 1*circshift(eye(SIZE),-1);

%initial amount of each species
x0 = zeros(SIZE,1);
x0(1) = 1;

%define reaction rates
a = 1.5;

%interpolation time vector and Tmax
Delta = 1/5000;
N = 2^15;
t_query = 0:Delta:(N-1)*Delta;
Tmax = t_query(end) + 10000*Delta;

%number of trials
bigN = 1;

%frequency vector
kkk = (-N/2:N/2-1);
fk = 1/(N*Delta)*kkk*2*pi;

%number of repititions for the power spectra
M = 200;

%yes
power_x = 0;
power_x_Q = 0;


%form the matrix
A = (-eye(SIZE)*a + a*circshift(eye(SIZE),-1))';
[v,lambda] = eig(A);
[w,lambda] = eig(A.');
lambda = diag(lambda)
imag(lambda)./real(lambda)
W = w(:,1)
P0 = v(:,3)/sum(v(:,3))


%loop for power spectra
for mmm = 1:M
    
    %solution vector
    SOLBOI = zeros(SIZE,length(t_query));
    SOLQ = zeros(1,length(t_query));
    
    %loopin it
    for i=1:bigN
        
        %define propensity function
        prop = @(x) a*[x(1);x(2);x(3)];
        
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
            alpha=prop(X(:,count));
            xi=rand;
            s=-1/sum(alpha)*log(1-xi);
            
            %draw the reaction
            phi=0;
            ell=0;
            eta=rand;
            while (phi < eta)
                ell=ell+1;
                phi=phi + alpha(ell)/sum(alpha);
            end
            
            %update
            t=t+s;
            X(:,count+1)=X(:,count) + S(:,ell);
            tspan(count+1)=t;
            
            %update count
            count = count + 1;
        end
        
        

        
        %do the smush smush
        X_sol = zeros(1,length(X));
        for pp=1:length(X_sol)
            X_sol(pp) = find(X(:,pp)==1);
        end
        
        %interpolate the data onto an evenly spaced time grid
        
        X_interped = interp1(tspan,X_sol,t_query,'previous');
        
        
        
        %do it
        for pp = 1:length(t_query)
            idx = X_interped(pp);
            SOLBOI(idx, pp) = SOLBOI(idx, pp) + 1;
        end
        
    end
    
    %normalize the solution to reflect the percentage of realizations in a
    %state
    SOLBOI = SOLBOI/bigN;
    
    
    %do the Q thing
    for yy=1:length(t_query)
        SOLQ(yy) = dot(W,SOLBOI(:,yy));
        SOLQ(yy) = SOLQ(yy)/(sum(abs(SOLQ(yy))));
    end
    
    %fft
    Fc_x = fftshift(fft(SOLBOI(2,:) - 0*mean(SOLBOI(2,:))));
    Fc_Qx = fftshift(fft(SOLQ - 0*mean(SOLQ)));
    
    %power spectra
    power_x = power_x + abs(Fc_x).^2;
    power_x_Q = power_x_Q + abs(Fc_Qx).^2;
end

%normalize those rascals
power_x = power_x/M;
power_x_Q = power_x_Q/M;

%exact
%exact = -2*real(lambda)./(real(lambda)^2 + (fk-imag(lambda)).^2);
%maxmax = max(exact);

%normalize
power_x_Q = power_x_Q/N*Delta;
%pdaddy_x_Q = pdaddy_x_Q/max(pdaddy_x_Q)*maxmax;

power_x = power_x/N*Delta;
%pdaddy_x = pdaddy_x/max(pdaddy_x)*maxmax;
power_x(power_x>0.4) = nan;

%plot power spectra
figure(7346)
hold on
semilogy(fk,(power_x),'-','color',[255, 182, 193 256]/256,'linewidth',2)
xlabel('frequency \nu')
ylabel('S_1(\nu)')
box on
axis square
set(gca,'FontSize',15)
xlim([-150 150])
ylim([10^-8 10^2])
set(gca, 'YScale', 'log')





