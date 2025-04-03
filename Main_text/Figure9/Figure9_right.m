%% setup

%mr clean
clc
clf

%parameters
a = 2;

%parameters
w2 = 1.9:.001:2.1;
gamma = 0:.001:.1;

%keep track of where the bifurcation appears for each (b,gamma) pair
dif = zeros(1,length(gamma));

%keep track of the TONGUE
TONGUE = zeros(length(gamma),length(w2));

%size of mah boi
n = 3;


%% form the matrix

%loopin' it
for qq = 1:length(w2)
    for pp = 1:length(gamma)
        
        %well
        GAMMA = gamma(pp);
        bb = w2(qq);
        
        %joint matrix without coupling
        A = -eye(n)*a + a*circshift(eye(n),-1);
        B = -eye(n)*bb + bb*circshift(eye(n),-1);
        well = kron(B,eye(n)) + kron(eye(n),A);
        
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
        well = well + BETA + ALPHA;
        
        
        %% diagonalize
        
        %diagonalize
        [~,lambda] = eig(well);
        lambda = diag(lambda);
        lambda = lambda(imag(lambda)>0);
        [~,I] = sort(real(lambda),'descend');
        lambda = lambda(I);
        lambda = lambda(1:2);
        dif(pp) = abs(lambda(1)-lambda(2));
    end
    
    %find min of dif to figure out where bifurcation occurs
    [~,I] = min(dif);
    
    %update the tongue
    TONGUE(I:end,qq) = 1;
    
end

%% visualize

%Define custom colormap
custom_colormap = [
    0 0 0;        
    [255, 182, 193]/256      
    ];

%plot the tongue
figure(1)
hold on
imagesc(w2-a,gamma,TONGUE)
plot(w2-a,sqrt(3)/2*abs(w2-a),'-','color',[1 0.6 0],'linewidth',4)
set(gca,'Ydir','normal')
xlabel('\tau')
ylabel('\kappa')
set(gca,'fontsize',15)
box on
axis square
colormap(custom_colormap)
xlim([min(w2-a) max(w2-a)])
ylim([min(gamma) max(gamma)])




