function Wave
c = 1; % constant in the wave equation
k = 2; % constant in exact solution
U = @ (x,t) cos(k*pi*t).*cos(k*pi*x); % wave equation for odd k

N = 4; % polynomial order
K = 12; % B-spline sections

% knot sequence
VX = linspace(-1,1,K+1); % internal knots on domain
at = -1*ones(1,N);
ct = 1*ones(1,N);
t = [at VX ct];
dt = [VX(1)*ones(1,N) VX(2:end-1) VX(end)*ones(1,N)];

% mapping
h = @(r) repmat(diff(VX),length(r),1);
map = @(r) reshape(h(r).*(repmat(r,1,K)+1)/2 + repmat(VX(1:end-1),length(r),1),length(r)*K,1);

% interpolation points over each interval for initial basis
[rq, wq] = JacobiGQ(0,0,N);
rBq = map(rq);
wBq = repmat(wq,1,K).*h(rq)/2;
wBq = wBq(:);

% interpolation points for final basis
r = zeros(n+k,1);
for i = 1:n+k
    r(i) = mean(t((i+1):(i+n))); % greville abscissae
end

DBr = DBrClosed(n,k); % spline coefficients to spline derivative coefficients
Dt = 1/(200*n*k); % time step
A = 1; % derivative constant
test = 1; % 1 for Forward Euler. 2 for ode45

%Fekete points
while A > 10^-12
    Br = splineDerv(VX,r,t,n,DBr); % dphi_i/dr_i
    Br(1) = 0; Br(end) = 0; % sets boundary conditions
    [r,A] = timeStep(Br,r,Dt,test); % adjusts r
end

% replacing initial basis with final basis
qBq = basisTransition(VX,N,rBq);
qBf = basisTransition(VX,N,rf);

% evaluate splines at interpolation points
Bq = bspline_basismatrix(N+1,t,rf); % Bq = Bq*inv(B)
dB = bspline_basismatrix(N,dt,rf);
DBr = DBrClosed(N,K);
DBq = dB*DBr; % DBr: spline coefficients to spline derivative coefficients

% transform to fekete nodal basis
wBf = qBf'\(qBq'*wBq); % desired results

dt = 1/89600; % time step

% compute matrices
M = Bq'*diag(wBf)*Bq; % non-diagonal mass matrix
ML = diag(wBf); % diagonal mass matrix
KK = DBq'*diag(wBf)*DBq; % stiffness matrix

%Apply corrector method
r = 80;
A = M/ML - eye(size(M));
add = 0;
for i = 0:r-1
    add = add + (-A)^i;
end

% initial conditions
u0 = @(x) cos(k*pi*x);
a = Bq\u0(rf); %(Bq'*diag(wBq)*Bq)\(Bq'*diag(wBq)*(u0(rBq)));
a = repmat(a,[1,2]);

% Animates collocation results through time
xx = linspace(-1,1,500)';
B = bspline_basismatrix(N+1,t,xx);
check = 896;
count = 1; % used for movie
for i = 0:dt:1
    if check == 896
        check = 0;
        u = B*a(:,2);
        plot(xx,U(xx,i));
        hold on
        plot(xx,u)
        ylim([-2 2]) % so movie frames don't auto-adjust. May need to change depending on solution range
        hold off
        F(count) = getframe;
        count = count + 1;
    end
    check = check + 1;

    % FE time leap-frogging scheme
    a(:,1) = (-ML\add*KK*a(:,2))*(c*dt)^2 + 2*a(:,2) - a(:,1); % explicit time leapfrogging method
    a = circshift(a, [0,1]);
end
movie(F)

% calculates L^2 norm
u = Bq*a(:,2);
error = sqrt(sum((U(rf,i) - u).^2.*wBf))
return