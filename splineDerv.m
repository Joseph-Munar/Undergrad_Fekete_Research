function u = splineDerv(VX,r,t,N,DBr)
% represent derivatives of splines as a spline of lower degree and
% continuity (lower continuity = doubled interior knots)
dt = [VX(1)*ones(1,N) VX(2:end-1) VX(end)*ones(1,N)];

%Pulls out dphi_i/dr_i for gradient ascent algorithm
B = bspline_basismatrix(N+1,t,r);
dB = bspline_basismatrix(N,dt,r);
u = dB*DBr/B; % derivative matrix
u = diag(u);
return