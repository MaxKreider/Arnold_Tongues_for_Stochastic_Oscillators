%% setup

%well
format long

%mr clean
clc

%starting spot
M = 0;

%choose diagonal row
N = -1;

%cutoff point to approximate infinite continued fraction as a finite one
cutoff = 200;

%fixed system parameters
D = 0.3;                            %change noise as needed
w1 = 2;

%variable system parameters   
w2 = 1.999:.00001:2.001;            %for figure 10 (right)
gamma = 0:.00001:.001;

%variable system parameters
w2 = 1.9:.001:2.1;                  %for figure 10 (left and middle)
gamma = 0:.001:.1;

%solution matrix
TONGUE = zeros(length(gamma),length(w2));


%% explicit matrix method

%loop over all parameter space
for jj = 1:length(w2)
    jj
    for ii = 1:length(gamma)
                
        %initialize matrix
        my_matrix_friend = zeros(2*cutoff+1,2*cutoff+1);
        
        %form the matrix
        count = 1;
        for j=-cutoff+M:cutoff+M
            
            %diagonal entries
            my_matrix_friend(count,count) = Q_mid(N,j,D,w1,w2(jj));
            
            %super/subdiagonal entries
            my_matrix_friend(count,count+1) = Q_plus(N,j,gamma(ii));
            my_matrix_friend(count+1,count) = Q_minus(N,j+1,gamma(ii));
            
            %update count
            count = count+1;
        end
        
        if mod(N,2) == 0
            my_matrix_friend = my_matrix_friend(1:2*cutoff+1,1:2*cutoff+1);
        else
            my_matrix_friend(end,end) = Q_mid(N,cutoff+M+1,D,w1,w2(jj));
        end
        
        %diagonalize
        lambda_boi = eig(my_matrix_friend);
        
        %sort
        [~,temp] = sort(-real(lambda_boi));
        lambda_boi = lambda_boi(temp);
        
        %get the eigenvalues to compare
        lambda_1 = lambda_boi(1);
        lambda_2 = lambda_boi(2);
        
        %compare
        if abs(imag(lambda_1) - imag(lambda_2)) < 1e-6
            TONGUE(ii:end,jj) = 1;
            break
        end
        
    end
end


%% visualize

% Define custom colormap
custom_colormap = [
    0 0 0;         
    0.6 0.6 0.6    
];

%plot
figure(1)
hold on
imagesc(w2-w1,gamma,TONGUE)
plot(w2-w1,abs((w2-w1)),'-','color',[1 0.6 0],'linewidth',4)
plot(w2-w1,abs((w2-w1)/2),'c--','linewidth',2)
set(gca, 'YDir','normal')
xlabel('\tau')
ylabel('\kappa')
set(gca,'fontsize',15)
box on
axis square
xlim([w2(1)-w1, w2(end)-w1])
ylim([0 gamma(end)])
colormap(custom_colormap)


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
