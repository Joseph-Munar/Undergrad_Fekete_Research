function r = feketePoints(N,K,t,dt)
% first if-statement: generate Fekete points
% second if-statement: pulls Fekete points from Excel docs on Github page

if 1 % 1 to run. 0 to not run
% calculates Fekete points for B-splines of order N, knot vector t, and K
% internal knots

r = zeros(N+K,1);
for i = 1:N+K
    r(i) = mean(t((i+1):(i+N))); % greville abscissae
end

DBr = DBrClosed(N,K); % spline coefficients to spline derivative coefficients
Dt = 1/(200*N*K); % time step
A = 1; % derivative constant
test = 2; % 1 for Forward Euler. 2 for ode45

% Gradient Ascenscion Algorithm
while A > 10^-12
    Br = splineDerv(r,t,dt,N,DBr); % dphi_i/dr_i
    Br(1) = 0; Br(end) = 0; % sets boundary conditions
    [r,A] = timeStep(Br,r,Dt,test); % adjusts r
end
end

if 0
    X = [num2str(N),'KnotFeketeMeshes.xlsx'];
    R = xlsread(X);
    r = R(K,1:(K+N+3))';
end
return
