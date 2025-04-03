%% setup

%mr clean
clc
clf

%parameters
w1 = 2;
w2 = 2.5;
gamma = 0.09;

%time
dt = 0.001;
Tmax = 5000;
t = 0:dt:Tmax;

%solution vectors
X = zeros(2,length(t));

%initial conditions
X(:,1) = rand(2,1)*2*pi;


%% do the computations for the histogram

%euler 
for pp=1:length(t)-1
    
    %integrate
    X(:,pp+1) = X(:,pp) + dt*[w1 + gamma*sin(X(2,pp)-X(1,pp)); w2 + gamma*sin(X(1,pp)-X(2,pp))];
    
    %enfornce periodicity
    X(:,pp+1) = mod(X(:,pp+1),2*pi);
end


%% chop off the transients

%plot chop
plot_chop = 1:33000;

%chop those rascals
X_plot = X(:,plot_chop);
X_hist = X(:,10000:end);


%% make a plot

% Use tiled layout with 2 rows, 3 columns
tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

nexttile(4)
hold on
scatter(t(plot_chop), X_plot(1,:), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
scatter(t(plot_chop), X_plot(2,:), 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7], 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
xlabel('time t')
ylabel('x(t)')
set(gca,'fontsize',15)
xlim([0 33])
ylim([0 2*pi])


%% make a histogram

%yeah
nexttile(5)
histo=histogram2(X_hist(1,:),X_hist(2,:),[100 100],'Normalization','pdf','DisplayStyle','tile','ShowEmptyBins','on','BinWidth',[2*pi/100 2*pi/100]);
shading interp
ax1 = gca; 
colormap(ax1, hot) 
xlabel('x')
ylabel('y')
axis square
set(gca,'FontSize',15)
xlim([0 2*pi])
ylim([0 2*pi])
colorbar
caxis([0 max(histo.Values(:))])


%% setup

%parameters
w1 = 2;
w2 = 2.5;
gamma = 0.26;

%time
dt = 0.001;
Tmax = 50000;
t = 0:dt:Tmax;

%solution vectors
X = zeros(2,length(t));

%initial conditions
X(:,1) = rand(2,1)*2*pi;


%% do the computations for the histogram

%euler <3
for pp=1:length(t)-1
    
    %integrate
    X(:,pp+1) = X(:,pp) + dt*[w1 + gamma*sin(X(2,pp)-X(1,pp)); w2 + gamma*sin(X(1,pp)-X(2,pp))];
    
    %enfornce periodicity
    X(:,pp+1) = mod(X(:,pp+1),2*pi);
end


%% chop off the transients

%plot chop
plot_chop = 1:33000;

%chop those rascals
X_plot = X(:,plot_chop);
X_hist = X(:,10000:end);


%% make a plot

nexttile(1)
hold on
scatter(t(plot_chop), X_plot(1,:), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
scatter(t(plot_chop), X_plot(2,:), 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7], 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
xlabel('time t')
ylabel('x(t)')
set(gca,'fontsize',15)
xlim([0 33])
ylim([0 2*pi])


%% make a histogram

%yeah
nexttile(2)
histo=histogram2(X_hist(1,:),X_hist(2,:),[100 100],'Normalization','pdf','DisplayStyle','tile','ShowEmptyBins','on','BinWidth',[2*pi/100 2*pi/100]);
shading interp
ax2 = gca; 
colormap(ax2, hot) 
xlabel('x')
ylabel('y')
axis square
set(gca,'FontSize',15)
xlim([0 2*pi])
ylim([0 2*pi])
caxis([0 max(histo.Values(:))])
cb = colorbar;


%% tongue

%parameters
w2 = 1.9:.001:2.1;
gamma = 0:.001:.1;

%solution vector
TONGUE = zeros(length(gamma),length(w2));


%% analyze eigenvalue bifurcations

%loopin it
for jj=1:length(w2)
    for ii=1:length(gamma)
        
        %check for equality of the imaginary parts
        if (gamma(ii) > abs(w2(jj)-2)/2)
            TONGUE(ii:end,jj) = 1;
            break
        end
        
    end
end


%% visualize

%colormap
custom_colormap = [0 0 0; 0.3 0 0.5];

%plot
nexttile(3, [2,1])
hold on
imagesc(w2-w1,gamma,TONGUE)
%plot(w2-w1,abs((w2-w1)/2),'-','color',[1 0.6 0],'linewidth',4)
set(gca, 'YDir','normal')
xlabel('\tau')
ylabel('\kappa')
set(gca,'fontsize',12)
box on
axis square
xlim([-.1 .1])
ylim([0 gamma(end)])
ax3 = gca; 
colormap(ax3, custom_colormap) 

set(gcf,'position',[181.8,275.4,1125.6,420.0000000000001])







