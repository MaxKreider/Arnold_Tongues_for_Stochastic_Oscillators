%% setup

%mr clean
clc

%parameters
D = .1;
k = .1;
w1 = 2;
w2 = 2.5;       %w2 = omega + tau
kappa = .1;    

%construct A
A = [-k-kappa, w1, kappa, 0;
    -w1, -k-kappa, 0, kappa;
    kappa, 0, -k-kappa, w2;
    0, kappa, -w2, -k-kappa];

%construct B
BB = 2*D*eye(4);

%covariance matrix
Sigma = lyap(A,BB);

%inverse of covariance matrix
Xi = inv(Sigma);

%setup for psd do-da
Delta = 1/500;
N = 2^17;
t = 0:Delta:(N-1)*Delta;

%frequency vector
kkk = (-N/2:N/2-1);
fk = 1/(N*Delta)*kkk*2*pi;

%number of trials
M = 4000;

%yes
power_x = 0;
power_y = 0;
power_cross = 0;
dad_cross = 0;
dad_x = 0;
dad_y = 0;


%% find the Q functions

%find left eigenfunctions and eigenvalues
[W,lambda] = eig(A.');
Qx = W(:,1);
Qy = W(:,3);
LAMBDA = [lambda(1,1);lambda(3,3)];


%% normalize the magnitude of the Q functions

%create a function
f1 = @(x1,y1,x2,y2) 1/(4*pi^2)*1/sqrt(det(Sigma))*(Qx(1)*x1+Qx(2)*y1+Qx(3)*x2+Qx(4)*y2).*conj(Qx(1)*x1+Qx(2)*y1+Qx(3)*x2+Qx(4)*y2).*exp(-1/2*(x1.*(Xi(4,1)*y2 + Xi(3,1)*x2 + Xi(1,1)*x1))).*exp(-1/2*(x2.*(Xi(2,3)*y1 + Xi(3,3)*x2 + Xi(1,3)*x1))).*exp(-1/2*(y1.*(Xi(4,2)*y2 + Xi(2,2)*y1 + Xi(3,2)*x2))).*exp(-1/2*(y2.*(Xi(4,4)*y2 + Xi(2,4)*y1 + Xi(1,4)*x1)));
f2 = @(x1,y1,x2,y2) 1/(4*pi^2)*1/sqrt(det(Sigma))*(Qy(1)*x1+Qy(2)*y1+Qy(3)*x2+Qy(4)*y2).*conj(Qy(1)*x1+Qy(2)*y1+Qy(3)*x2+Qy(4)*y2).*exp(-1/2*(x1.*(Xi(4,1)*y2 + Xi(3,1)*x2 + Xi(1,1)*x1))).*exp(-1/2*(x2.*(Xi(2,3)*y1 + Xi(3,3)*x2 + Xi(1,3)*x1))).*exp(-1/2*(y1.*(Xi(4,2)*y2 + Xi(2,2)*y1 + Xi(3,2)*x2))).*exp(-1/2*(y2.*(Xi(4,4)*y2 + Xi(2,4)*y1 + Xi(1,4)*x1)));

%integrate
I1 = integralN(f1,-10,10,-10,10,-10,10,-10,10,'AbsTol',1e-5,'RelTol',1e-3);
I2 = integralN(f2,-10,10,-10,10,-10,10,-10,10,'AbsTol',1e-5,'RelTol',1e-3);

%normalize
Qx = Qx/sqrt(I1);
Qy = Qy/sqrt(I2);


%% gauge transformation

%create a function
f1a = @(x1,y1,x2,y2) 1/(4*pi^2)*1/sqrt(det(Sigma))*(Qx(1)*x1+Qx(2)*y1+Qx(3)*x2+Qx(4)*y2).*conj(Qy(1)*x1+Qy(2)*y1+Qy(3)*x2+Qy(4)*y2).*exp(-1/2*(x1.*(Xi(4,1)*y2 + Xi(3,1)*x2 + Xi(1,1)*x1))).*exp(-1/2*(x2.*(Xi(2,3)*y1 + Xi(3,3)*x2 + Xi(1,3)*x1))).*exp(-1/2*(y1.*(Xi(4,2)*y2 + Xi(2,2)*y1 + Xi(3,2)*x2))).*exp(-1/2*(y2.*(Xi(4,4)*y2 + Xi(2,4)*y1 + Xi(1,4)*x1)));
f2a = @(x1,y1,x2,y2) 1/(4*pi^2)*1/sqrt(det(Sigma))*(Qy(1)*x1+Qy(2)*y1+Qy(3)*x2+Qy(4)*y2).*conj(Qx(1)*x1+Qx(2)*y1+Qx(3)*x2+Qx(4)*y2).*exp(-1/2*(x1.*(Xi(4,1)*y2 + Xi(3,1)*x2 + Xi(1,1)*x1))).*exp(-1/2*(x2.*(Xi(2,3)*y1 + Xi(3,3)*x2 + Xi(1,3)*x1))).*exp(-1/2*(y1.*(Xi(4,2)*y2 + Xi(2,2)*y1 + Xi(3,2)*x2))).*exp(-1/2*(y2.*(Xi(4,4)*y2 + Xi(2,4)*y1 + Xi(1,4)*x1)));

%integrate
I1a = integralN(f1a,-10,10,-10,10,-10,10,-10,10,'AbsTol',1e-5,'RelTol',1e-3);
I2a = integralN(f2a,-10,10,-10,10,-10,10,-10,10,'AbsTol',1e-5,'RelTol',1e-3);

%rotation normalization
factor = log(I1a/I2a);
if imag(factor)<0
    factor = factor + 2*pi*sqrt(-1);
end
alpha = -sqrt(-1)/2*factor;

%normalize
Qy = Qy*exp(sqrt(-1)*alpha);


%% compute the bracket term

%create a function
f = @(x1,y1,x2,y2) 1/(4*pi^2)*1/sqrt(det(Sigma))*(Qx(1)*x1+Qx(2)*y1+Qx(3)*x2+Qx(4)*y2).*conj(Qy(1)*x1+Qy(2)*y1+Qy(3)*x2+Qy(4)*y2).*exp(-1/2*(x1.*(Xi(4,1)*y2 + Xi(3,1)*x2 + Xi(1,1)*x1))).*exp(-1/2*(x2.*(Xi(2,3)*y1 + Xi(3,3)*x2 + Xi(1,3)*x1))).*exp(-1/2*(y1.*(Xi(4,2)*y2 + Xi(2,2)*y1 + Xi(3,2)*x2))).*exp(-1/2*(y2.*(Xi(4,4)*y2 + Xi(2,4)*y1 + Xi(1,4)*x1)));

%integrate
result = integralN(f,-10,10,-10,10,-10,10,-10,10,'AbsTol',1e-5,'RelTol',1e-3)


%% compute numerical trajectories

parfor qq=1:M
    
    qq
    
    %solution vectors
    Qxt = zeros(1,length(t));
    Qyt = zeros(1,length(t));
    XXt = zeros(1,length(t));
    YYt = zeros(1,length(t));
    
    %trajectory initial value (uniformly random on the torus)
    SOL = randn(4,1);
    
    %initial perturbed phase
    Qxt(1) = Qx(1)*SOL(1)+Qx(2)*SOL(2)+Qx(3)*SOL(3)+Qx(4)*SOL(4);
    Qyt(1) = Qy(1)*SOL(1)+Qy(2)*SOL(2)+Qy(3)*SOL(3)+Qy(4)*SOL(4);
    XXt(1) = SOL(1);
    YYt(1) = SOL(3);
    
    %compute phases
    for pp=1:length(t)-1
        
        %EM method
        SOL = SOL + Delta*A*SOL + sqrt(Delta)*sqrt(2*D)*randn(4,1);
        
        %compute Q function at that point
        Qxt(pp+1) = Qx(1)*SOL(1)+Qx(2)*SOL(2)+Qx(3)*SOL(3)+Qx(4)*SOL(4);
        Qyt(pp+1) = Qy(1)*SOL(1)+Qy(2)*SOL(2)+Qy(3)*SOL(3)+Qy(4)*SOL(4);
        XXt(pp+1) = SOL(1);
        YYt(pp+1) = SOL(3);
        
    end
    
    %take fft
    Fc_x = fftshift(fft(Qxt));
    Fc_y = fftshift(fft(Qyt));
    Fc_X = fftshift(fft(XXt));
    Fc_Y = fftshift(fft(YYt));
    
    %power spectra
    power_x = power_x + abs(Fc_x).^2;
    power_y = power_y + abs(Fc_y).^2;
    dad_x = dad_x + abs(Fc_X).^2;
    dad_y = dad_y + abs(Fc_Y).^2;
    
    %cross spectra
    power_cross = power_cross  + (Fc_x).*conj(Fc_y);
    dad_cross = dad_cross  + (Fc_X).*conj(Fc_Y);
end

%normalize
power_x = power_x/M;
power_y = power_y/M;
power_cross = power_cross/M;
dad_cross = dad_cross/M;
dad_x = dad_x/M;
dad_y = dad_y/M;

%exact solutions for power spectra
mu_x = real(LAMBDA(1));
omega_x = imag(LAMBDA(1));
mu_y = real(LAMBDA(2));
omega_y = imag(LAMBDA(2));
exact_x = 2*abs(mu_x)./(mu_x^2+(fk-omega_x).^2);
exact_y = 2*abs(mu_y)./(mu_y^2+(fk-omega_y).^2);

%exact solutions for cross spectra
expression = 1./((LAMBDA(1))-sqrt(-1)*fk) + 1./(conj(LAMBDA(2))+sqrt(-1)*fk);
expression = expression*(-result);


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
ylim([0 12])
xlabel('frequency \nu')
ylabel('S_1(\nu)')
box on
axis square
set(gca,'fontsize',15)

%plot cross spectra
figure(2)
hold on
plot(fk,real(dad_cross)/N*Delta,'-','color',[.6 .6 .6 .6],'linewidth',10)
plot(fk,-real(power_cross)/N*Delta,'-','color','#006400','linewidth',10)
plot(fk,-real(expression),'-','color',[0 1 0],'linewidth',3)
xlim([1 3.5])
ylim([0 6.5])
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
ylim([-1.5 1.5])
xlabel('frequency \nu')
ylabel('Im($S_{\lambda_x,\lambda_y}$)', 'Interpreter', 'latex');
box on
axis square
set(gca,'fontsize',15)


%% sigh

function q = integralN(f,varargin)
%Iterate INTEGRAL2 and INTEGRAL3 functions to calculate integrals of order
%4, 5, and 6. For convenience, the function is a straight pass-through to
%INTEGRAL, INTEGRAL2, and INTEGRAL3 for integrals of order 1, 2, and 3,
%respectively.
%
% Notes:
% *  This code has not been thoroughly tested.
% *  As with INTEGRAL, INTEGRAL2, and INTEGRAL3, your integrand MUST be
%    capable of accepting arrays and returning an array of the same size
%    with the value of the integrand computed "element-wise". In other
%    words, the method performs function evaluations in "batches". If your
%    integrand function only works with scalar arguments, instead of
%    passing f to this function directly, pass in...
%    ...for a 4D integral: @(x,y,z,w)arrayfun(f,x,y,z,w)
%    ...for a 5D integral: @(x,y,z,w,v)arrayfun(f,x,y,z,w,v)
%    ...for a 6D integral: @(x,y,z,w,v,u)arrayfun(f,x,y,z,w,v,u)
%    The same principle applies to functions that are used to define the
%    region of integration.
% *  Iterated adaptive quadrature is not a very efficient way of
%    approaching high order integrals. If speed is important, consider
%    sparse grid and Monte Carlo methods, among others, instead. Expect
%    this function to be slow on order 4 integrals, painfully slow on order
%    5, and excruciatingly slow on order 6. When using infinite limits (or
%    the 'iterated' method to deal with non-smoothness of the integrand) it
%    will take even longer. For example, the 4D hypersphere problem below
%    took less than 2 seconds on my machine, the 5D hypersphere problem
%    took over 1.5 minutes, and the 6D hypersphere problem took over 24
%    minutes. The example integrating over all of R^4 took nearly 8
%    minutes.
% *  I recommend that you loosen the tolerances considerably, since for
%    technical reasons this method is likely to produce more accuracy than
%    requested on smooth problems. See the examples below. Always supply
%    both RelTol and AbsTol, never just one or the other. For an
%    explanation of these tolerances see
%
%    http://blogs.mathworks.com/loren/2013/12/26/double-integration-in-matlab-understanding-tolerances/
%
% Examples:
%    % Calculate the "volume" of the unit hypersphere.
%    xmin = -1;
%    xmax = 1;
%    ymin = @(x)-sqrt(1 - x.*x);
%    ymax = @(x) sqrt(1 - x.*x);
%    zmin = @(x,y)-sqrt(1 - x.*x - y.*y);
%    zmax = @(x,y) sqrt(1 - x.*x - y.*y);
%    wmin = @(x,y,z)-sqrt(1 - x.*x - y.*y - z.*z);
%    wmax = @(x,y,z) sqrt(1 - x.*x - y.*y - z.*z);
%
%    % 4D
%    f4 = @(x,y,z,w)ones(size(x));
%    q4 = integralN(f4,xmin,xmax,ymin,ymax,zmin,zmax,wmin,wmax,'AbsTol',1e-5,'RelTol',1e-3)
%    error = q4 - pi^2/2
%
%    % 5D
%    vmin = @(x,y,z,w)-sqrt(1 - x.*x - y.*y - z.*z  - w.*w);
%    vmax = @(x,y,z,w) sqrt(1 - x.*x - y.*y - z.*z  - w.*w);
%    f5 = @(x,y,z,w,v)ones(size(x));
%    q5 = integralN(f5,xmin,xmax,ymin,ymax,zmin,zmax,wmin,wmax,vmin,vmax,'AbsTol',1e-5,'RelTol',1e-3)
%    error = q5 - 2*pi^(5/2)/(gamma(5/2)*5)
%
%    % 6D
%    umin = @(x,y,z,w,v)-sqrt(1 - x.*x - y.*y - z.*z  - w.*w - v.*v);
%    umax = @(x,y,z,w,v) sqrt(1 - x.*x - y.*y - z.*z  - w.*w - v.*v);
%    f6 = @(x,y,z,w,v,u)ones(size(x));
%    q6 = integralN(f6,xmin,xmax,ymin,ymax,zmin,zmax,wmin,wmax,vmin,vmax,umin,umax,'AbsTol',1e-5,'RelTol',1e-3)
%    error = q6 - pi^3/6
%
%    % Calculate a rapidly decaying 4D function over R^4.
%    f = @(x,y,z,w)exp(-x.*x - y.*y - z.*z - w.*w);
%    q = integralN(f,-inf,inf,-inf,inf,-inf,inf,-inf,inf,'AbsTol',1e-3,'RelTol',1e-3)
%
% Author:  Michael Hosea
% Date:  2014-09-25
% Copyright 2014 The MathWorks, Inc.

narginchk(3,inf);
if ~isa(f,'function_handle') || nargin(f) == 0
    error('integralN:BadIntegrand', ...
        ['First argument must be a function handle accepting one or more inputs.\n', ...
        'To integrate a constant c, integrate @(x,y)c*ones(size(x)) for 2D,\n', ...
        '@(x,y,z)c*ones(size(x)) for 3D, etc.']);
end
switch(nargin(f))
    case 1
        q = integral(f,varargin{:});
    case 2
        q = integral2(f,varargin{:});
    case 3
        q = integral3(f,varargin{:});
    case 4
        q = integral4(f,varargin{:});
    case 5
        q = integral5(f,varargin{:});
    case 6
        q = integral6(f,varargin{:});
    otherwise
        error('integralN:NTooLarge','N >= 7 is not supported.');
end

%--------------------------------------------------------------------------
end

function q = integral4(f,varargin)
% integral4(f) = integral2(integral2(f)).
narginchk(9,inf);
zmin = varargin{5};
zmax = varargin{6};
wmin = varargin{7};
wmax = varargin{8};
anyisinf = false;
if ~isa(zmin,'function_handle')
    anyisinf = anyisinf || isinf(zmin(1));
    zmin = @(x,y)zmin(1)*ones(size(x));
end
if ~isa(zmax,'function_handle')
    anyisinf = anyisinf || isinf(zmax(1));
    zmax = @(x,y)zmax(1)*ones(size(x));
end
if ~isa(wmin,'function_handle')
    anyisinf = anyisinf || isinf(wmin(1));
    wmin = @(x,y,z)wmin(1)*ones(size(x));
end
if ~isa(wmax,'function_handle')
    anyisinf = anyisinf || isinf(wmax(1));
    wmax = @(x,y,z)wmax(1)*ones(size(x));
end
if anyisinf
    method_override = {'method','iterated'};
else
    method_override = {};
end
inner = @(x,y)integral2( ...
    @(z,w)f(x*ones(size(z)),y*ones(size(z)),z,w), ...
    zmin(x,y), ...
    zmax(x,y), ...
    @(z)wmin(x*ones(size(z)),y*ones(size(z)),z), ...
    @(z)wmax(x*ones(size(z)),y*ones(size(z)),z), ...
    varargin{9:end},method_override{:});
q = integral2( ...
    @(xv,yv)arrayfun(inner,xv,yv), ...
    varargin{1:4},varargin{9:end});

%--------------------------------------------------------------------------
end

function q = integral5(f,varargin)
% integral5(4) = integral3(integral2(f)).
narginchk(11,inf);
wmin = varargin{7};
wmax = varargin{8};
vmin = varargin{9};
vmax = varargin{10};
anyisinf = false;
if ~isa(wmin,'function_handle')
    anyisinf = anyisinf || isinf(wmin(1));
    wmin = @(x,y,z)wmin(1)*ones(size(x));
end
if ~isa(wmax,'function_handle')
    anyisinf = anyisinf || isinf(wmax(1));
    wmax = @(x,y,z)wmax(1)*ones(size(x));
end
if ~isa(vmin,'function_handle')
    anyisinf = anyisinf || isinf(vmin(1));
    vmin = @(x,y,z,w)vmin(1)*ones(size(x));
end
if ~isa(vmax,'function_handle')
    anyisinf = anyisinf || isinf(vmax(1));
    vmax = @(x,y,z,w)vmax(1)*ones(size(x));
end
if anyisinf
    method_override = {'method','iterated'};
else
    method_override = {};
end
inner = @(x,y,z)integral2( ...
    @(w,v)f(x*ones(size(w)),y*ones(size(w)),z*ones(size(w)),w,v), ...
    wmin(x,y,z), ...
    wmax(x,y,z), ...
    @(w)vmin(x*ones(size(w)),y*ones(size(w)),z*ones(size(w)),w), ...
    @(w)vmax(x*ones(size(w)),y*ones(size(w)),z*ones(size(w)),w), ...
    varargin{11:end},method_override{:});
q = integral3( ...
    @(xv,yv,zv)arrayfun(inner,xv,yv,zv), ...
    varargin{1:6},varargin{11:end});

%--------------------------------------------------------------------------
end

function q = integral6(f,varargin)
% integral6(f) = integral4(integral2(f)).
narginchk(13,inf);
vmin = varargin{9};
vmax = varargin{10};
umin = varargin{11};
umax = varargin{12};
anyisinf = false;
if ~isa(vmin,'function_handle')
    anyisinf = anyisinf || isinf(vmin(1));
    vmin = @(x,y,z,w)vmin(1)*ones(size(x));
end
if ~isa(vmax,'function_handle')
    anyisinf = anyisinf || isinf(vmax(1));
    vmax = @(x,y,z,w)vmax(1)*ones(size(x));
end
if ~isa(umin,'function_handle')
    anyisinf = anyisinf || isinf(umin(1));
    umin = @(x,y,z,w,v)umin(1)*ones(size(x));
end
if ~isa(umax,'function_handle')
    anyisinf = anyisinf || isinf(umax(1));
    umax = @(x,y,z,w,v)umax(1)*ones(size(x));
end
if anyisinf
    method_override = {'method','iterated'};
else
    method_override = {};
end
inner = @(x,y,z,w)integral2( ...
    @(v,u)f(x*ones(size(v)),y*ones(size(v)),z*ones(size(v)),w*ones(size(v)),v,u), ...
    vmin(x,y,z,w), ...
    vmax(x,y,z,w), ...
    @(v)umin(x*ones(size(v)),y*ones(size(v)),z*ones(size(v)),w*ones(size(v)),v), ...
    @(v)umax(x*ones(size(v)),y*ones(size(v)),z*ones(size(v)),w*ones(size(v)),v), ...
    varargin{13:end},method_override{:});
q = integral4( ...
    @(xv,yv,zv,wv)arrayfun(inner,xv,yv,zv,wv), ...
    varargin{1:8},varargin{13:end});

%--------------------------------------------------------------------------
end