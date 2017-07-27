function splineFunction
% interpolates a continuous 1-D function using Fekete points and a Nth
% order B-spline basis with K knot intervals

%Insert function
f = @(x) (1-x).*cos(3*pi*(x+1)).*exp(-(1+x)); % insert function

N = 5; % spline degree
K = 12; %partitions of the spline (number of polynomial approximations needed)

% knot sequence
at = -1*ones(1,N);
ct = ones(1,N);
VX = linspace(-1,1,K+1); % internal knots on domain
t0 = [at VX ct]; % creates initial knot vector
t = t0;
re = linspace(-1,1,N+K)';
for i = 1:30
    Be = bspline_basismatrix(N+1, t, t0);
    t = (Be*re)'; % "smoothes' the knot out
end
dt = [VX(1)*ones(1,N) VX(2:end-1) VX(end)*ones(1,N)];

r = feketePoints(N,K,t,dt);
rp = -1:0.01:1; % finely spaced plotting grid

%Points on B-spline curve
B = bspline_basismatrix(N+1,t,r);
Bp = bspline_basismatrix(N+1,t,rp)/B;

%Plot control points and spline
figure
app = Bp*f(r);
hold on
plot(rp,app,'--','linewidth',2)
fplot(f, [-1,1]);
legend('Approximated Curve', 'Actual curve','Location','Best');
Title = ['n=', num2str(N), ' k=', num2str(K)];
title(Title)
hold off
error = abs(f(rp)'-app);
Lebesque = max(error);
figure
plot(rp,error)
Title = ['Error. Max = ', num2str(Lebesque)];
title(Title)
return