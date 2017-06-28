function splineFunction
%Insert function
f = @(x) (1-x).*cos(3*pi*(x+1)).*exp(-(1+x)); % insert function

n = 3; % spline degree
k = 12; %partitions of the spline (number of polynomial approximations needed)

% knot sequence
at = -1*ones(1,n);
ct = ones(1,n);
VX = linspace(-1,1,k+1); % internal knots on domain
t0 = [at VX ct]; % creates initial knot vector
t = t0;
re = linspace(-1,1,n+k)';
for i = 1:30
    Be = bspline_basismatrix(n+1, t, t0);
    t = (Be*re)'; % "smoothes' the knot out
end
VX = t(n+1:k+n+1); % interior points of smoothed knot

r = zeros(n+k,1);
for i = 1:n+k
    r(i) = mean(t((i+1):(i+n))); % greville abscissae
end

DBr = DBrClosed(n,k); % spline coefficients to spline derivative coefficients
dt = 1/(200*n*k); % time step
A = 1; % derivative constant
test = 1; % 1 for Forward Euler. 2 for ode45

%Fekete points
while A > 10^-12
    Br = splineDerv(VX,r,t,n,DBr); % dphi_i/dr_i
    Br(1) = 0; Br(end) = 0; % sets boundary conditions
    [r,A] = timeStep(Br,r,dt,test); % adjusts r
end
rp = -1:0.01:1; % finely spaced plotting grid

%Points on B-spline curve
B = bspline_basismatrix(n+1,t,r);
Bp = bspline_basismatrix(n+1,t,rp)/B;

%Plot control points and spline
figure
app = Bp*f(r);
hold on
plot(rp,app,'--','linewidth',2)
fplot(f, [-1,1]);
legend('Approximated Curve', 'Actual curve','Location','Best');
Title = ['n=', num2str(n), ' k=', num2str(k)];
title(Title)
hold off
error = abs(f(rp)'-app);
Lebesque = max(error);
figure
plot(rp,error)
Title = ['Error. Max = ', num2str(Lebesque)];
title(Title)
return