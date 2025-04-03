%% setup

%mr clean
clc
clf

%system paramters
k1 = .1;
k2 = .1;
w1 = 2;

%parameters
w2 = 1.9:.001:2.1;
gamma = 0:.001:.1;

%solution vector
TONGUE = zeros(length(gamma),length(w2));


%% create system

%coupling matrix
B = [-1 0  1  0;
    0 -1  0  1;
    1  0 -1  0;
    0  1  0 -1];


%% analyze eigenvalue bifurcations

%loopin it
for jj=1:length(w2)
    jj
    for ii=1:length(gamma)
        
        %uncoupled system
        A = [-k1 w1  0      0;
            -w1 -k1  0      0;
            0    0  -k2     w2(jj);
            0    0  -w2(jj) -k2];
        
        %full system
        L = A+gamma(ii)*B;
                
        %look at eigenvalues / eigenvectors
        [W,D] = eig(L.');
        
        %check for equality of the imaginary parts
        if abs(imag(D(1,1)) - imag(D(3,3))) < 1e-6
            TONGUE(ii:end,jj) = 1;
            break
        end
               
    end
end




%% visualize

% Define custom colormap
custom_colormap = [
    0 0 0;        
    0 1 0      
    ];

%plot
figure(1)
hold on
imagesc(w2-w1,gamma,TONGUE)
plot(w2-w1,abs((w2-w1)/2),'-','color',[1 0.6 0],'linewidth',4)
set(gca, 'YDir','normal')
xlabel('\tau')
ylabel('\kappa')
set(gca,'fontsize',15)
box on
axis square
xlim([-.1 .1])
ylim([0 gamma(end)])
colormap(custom_colormap)







