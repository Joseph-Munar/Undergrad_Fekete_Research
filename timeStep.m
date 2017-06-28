function [r,a] = timeStep(Vr,r, test)
%Gradient Ascension algorithm
dt = 0.001; % time step per iteration
if test == 1
    %Foward Euler approximation
    R = r;
    R = R + Vr*dt;
elseif test == 2
    %ode45 approximation
    tspan = [0 dt];
    [~, R] = ode45(@(r, phi) Vr, tspan, r);
    [m,~] = size(R);
    R = R(m,:)';
end
a = norm(R-r)/dt; % numerically checks derivative
r = R;
return