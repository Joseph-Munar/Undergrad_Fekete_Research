function feketePoints
f = @(x) 1./(1+10*x.^2); % insert funnction

N = 12; % degree of polynomial
r = linspace(-1,1,N+1)'; % equispaced nodes
a = 1; % derivative constraint
test = 1; % 1 for Forward Euler. 2 for ode45

while a > 10^-12 % run until derivative is very small
    V = Vandermonde1D(N,r);
    Vr = GradVandermonde1D(N,r)/V; % derivative matrix
    Vr = diag(Vr); %dphi_i/dr_i
    Vr(1) = 0; Vr(N+1) = 0; % sets boundaries
    [r,a] = timeStep(Vr,r,test);
end
V = Vandermonde1D(N,r);

rp = -1:0.01:1; % finely spaced plotting grid
Vp = Vandermonde1D(N,rp)/V;
figure
plot(rp,f(rp),'linewidth',2)
hold on
plot(rp,Vp*f(r),'--','linewidth',2)
plot(r, 0*r, 'O') % shows nodal points
legend('Exact function','Interpolant')
return