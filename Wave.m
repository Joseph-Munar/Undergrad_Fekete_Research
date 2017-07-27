function Wave
%Plots exact and approximated solution to the 1-D wave function using Nth
%order B-splines with K internal knots

c = 1; % constant in the wave equation
k = 2; % constant in exact solution
U = @ (x,t) cos(k*pi*t).*cos(k*pi*x);

N = 5; % polynomial order
K = 8; % B-spline sections

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
rf = feketePoints(N,K,t,dt);

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

Dt = 1/89600; % time step

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
a = Bq\u0(rf);
a = repmat(a,[1,2]);

% Animates collocation results through time
xx = linspace(-1,1,500)';
B = bspline_basismatrix(N+1,t,xx);
check = 896;
count = 1; % used for movie
for i = 0:Dt:1
    if check == 896
        check = 0;
        u = B*a(:,2);
        plot(xx,U(xx,i),'--');
        hold on
        plot(xx,u)
        ylim([-2 2]) % so movie frames don't auto-adjust. May need to change depending on solution range
        hold off
        F(count) = getframe;
        count = count + 1;
    end
    check = check + 1;

    % FE time leap-frogging scheme
    a(:,1) = (-ML\add*KK*a(:,2))*(c*Dt)^2 + 2*a(:,2) - a(:,1); % explicit time leapfrogging method
    a = circshift(a, [0,1]);
end
movie(F)

% calculates L^2 norm
u = Bq*a(:,2);
error = sqrt(sum((U(rf,i) - u).^2.*wBf))
return