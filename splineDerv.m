function u = splineDerv(r,t,dt,N,DBr)
%Pulls out dphi_i/dr_i for gradient ascent algorithm

B = bspline_basismatrix(N+1,t,r);
dB = bspline_basismatrix(N,dt,r);
u = dB*DBr/B; % derivative matrix
u = diag(u);
return