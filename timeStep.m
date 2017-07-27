function [r,a] = timeStep(Vr,r,Dt,test)
%Gradient Ascension Algorithm

if test == 1
    %Foward Euler approximation
    R = r;
    R = R + Vr*Dt;
elseif test == 2
    %ode45 approximation
    tspan = [0 Dt];
    [~, R] = ode45(@(r, phi) Vr, tspan, r);
    [m,~] = size(R);
    R = R(m,:)';
end
a = norm(R-r)/Dt; % numerically checks derivative
r = R;
return