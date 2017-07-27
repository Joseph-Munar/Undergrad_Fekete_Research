function qB = basisTransition(VX,N,rB)
% creates matrix qB of the form qB_i,j = integral(phi_i(x)*phi_j(x))dx

% knot sequence
VX = sort(repmat(VX(2:end-1), [1,N-1]));
at = -1*ones(1,2*N);
ct = ones(1,2*N);
t = [at VX ct];

qB = bspline_basismatrix(2*N,t,rB);
return