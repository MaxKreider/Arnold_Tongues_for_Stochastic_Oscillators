%% setup

%mr clean, gravy why you flow so mean
clc

%starting spot
M = 1;

%choose diagonal row
N = 1;

%cutoff point to approximate infinite continued fraction as a finite one
cutoff = 50;

%system parameters
D = .1;
w1 = 2;
w2 = 2.5;
gamma = 1;

%lambda
lambda =  -1.796113654986069 + 2.250000000000632i;

%solution vector
S = zeros(1,cutoff-1);

%indexing L :(
ccount = 1;

%spatial res
dt = 0.005;

%loop
for n=M:cutoff+M
    
    %% define numerator and denominator of CF in the style of Risken, pg. 227
    
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
    
    %% iterate
    
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
c = zeros(1,cutoff);
c(1) = 1;

%solve for coefficients
for j=1:cutoff
    c(j+1) = S(j)*c(j);
end

%% compute eigenfunction

%space vector
x = 0:dt:2*pi;
y = x';

%construct
Qf = 0;
count = 1;
for n=M:cutoff+M
    Qf = Qf + c(count).*exp(sqrt(-1)*(n*x+(N-n)*y));
    count = count+1;
end


%% backwards iteration

%reset
S = 0;

%indexing L :(
ccount = 1;

%loop
for n=M:-1:-cutoff+M
    
    %% define numerator and denominator of CF in the style of Risken, pg. 227
    
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
    
    %% iterate
    
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
c = zeros(1,cutoff);
c(1) = 1;

%solve for coefficients
for j=1:cutoff
    c(j+1) = S(j)*c(j);
end

%% compute eigenfunction

%space vector
x = 0:dt:2*pi;
y = x';

%construct
Qb = 0;
count = 2;
for n=M-1:-1:-cutoff+M
    Qb = Qb + c(count).*exp(sqrt(-1)*(n*x+(N-n)*y));
    count = count + 1;
end


%% save the data

%do it
Q_gamma = (Qf+Qb)*exp(sqrt(-1)*pi*0);
Qphase = angle(Q_gamma);
Qmag = abs(Q_gamma)/max(max(abs(Q_gamma)));


%% color map

%define a grid
npoints = 500;
[xxx, yyy] = meshgrid(linspace(-1, 1, npoints), linspace(-1, 1, npoints));

%compute polar coordinates
theta = atan2(yyy, xxx);          
rho   = sqrt(xxx.^2 + yyy.^2);   
mag   = min(rho, 1);          

%normalize phase from [-pi, pi] to [0,1] for hue
hue = (theta + pi) / (2*pi);

%saturation: vivid on the outside, washed out in the center
saturation = mag; 

%value (brightness): increases outward, but with a minimum baseline to prevent total darkness
value = 0.15 + 0.85 * mag; 

%compose the HSV image
hsv_img = cat(3, hue, saturation, value);

%convert HSV to RGB
rgb_img = hsv2rgb(hsv_img);

%put the `wheel' in colorwheel
mask = rho > 1;
rgb_img(repmat(mask, 1, 1, 3)) = 1;

%display the colorwheel
figure(1);
imshow(rgb_img, 'InitialMagnification', 'fit');
axis equal; 
set(gca, 'XColor', 'none', 'YColor', 'none');
colormap(hsv)
cbar = colorbar;
cbar.Ticks = linspace(0, 1, 5);
cbar.TickLabels = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'};
caxis([0,1]);
cbar.FontSize = 14;
cbar.Label.FontSize = 16;


%% do the mapping

%hue
hue = (Qphase + pi) / (2*pi);

%saturation
saturation = Qmag; 

%brightness
value = 0.15 + 0.85 * Qmag; 

%hsv image stuff
hsv_img = cat(3, hue, saturation, value);
rgb_img = hsv2rgb(hsv_img);

%display this rascally interpolated complex eigenfunction
figure(2);
imagesc(x,y,rgb_img); 
axis image;       
set(gca, 'YDir', 'normal');  
xlabel('x');
ylabel('y');
set(gca,'fontsize',15)


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
