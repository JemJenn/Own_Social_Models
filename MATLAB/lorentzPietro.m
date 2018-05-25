% lorenz.m
%
% Imlements the Lorenz System ODE
%   dx/dt = sigma*(y-x)
%   dy/dt = x*(rho-z)-y
%   dz/dt = xy-beta*z
%
% Inputs:
%   t - Time variable: not used here because our equation
%       is independent of time, or 'autonomous'.
%   x - Independent variable: has 3 dimensions
% Output:
%   dx - First derivative: the rate of change of the three dimension
%        values over time, given the values of each dimension
function dx = lorentzPietro(t, x)

% Standard constants for the Lorenz Attractor
alpha1 = 7;
alpha2 = 3;
beta1 = .3;
beta2 = .7;
b = 2;
c = 0.5;

% I like to initialize my arrays
dx = [0; 0; 0];

% The lorenz strange attractor
dx(1) = -x(1)*((b*x(1)*x(2) - b*beta2*x(3)*(x(1) - 1))*(x(2) - 1) - x(1)*((b*x(1)*x(2) - b*beta2*x(3)*(x(1) - 1))*(x(2) - 1) + x(2)*((x(1) - 1)*(beta1*c*(x(3) - 1) + beta1*beta2*x(3)*(b - c) - b*beta2*x(3)*(beta1 - 1) + beta1*c*x(3)*(beta2 - 1)) - x(1)*(c*(x(2) - 1) + alpha1*x(2)*(b - c)))) + x(2)*((x(1) - 1)*(beta1*c*(x(3) - 1) + beta1*beta2*x(3)*(b - c) - b*beta2*x(3)*(beta1 - 1) + beta1*c*x(3)*(beta2 - 1)) - x(1)*(c*(x(2) - 1) + alpha1*x(2)*(b - c))) - ((b*x(3)*(x(1) - 1) - b*beta1*x(1)*x(2))*(x(3) - 1) - x(3)*((c*(x(3) - 1) + alpha2*x(3)*(b - c))*(x(1) - 1) - x(1)*(beta2*c*(x(2) - 1) + beta1*beta2*x(2)*(b - c) - b*beta1*x(2)*(beta2 - 1) + beta2*c*x(2)*(beta1 - 1))))*(x(1) - 1));
dx(2) = x(2)*(x(1)*(c*(x(2) - 1) + alpha1*x(2)*(b - c)) - (x(1) - 1)*(beta1*c*(x(3) - 1) + beta1*beta2*x(3)*(b - c) - b*beta2*x(3)*(beta1 - 1) + beta1*c*x(3)*(beta2 - 1)) + (b*x(1)*x(2) - b*beta2*x(3)*(x(1) - 1))*(x(2) - 1) + x(2)*((x(1) - 1)*(beta1*c*(x(3) - 1) + beta1*beta2*x(3)*(b - c) - b*beta2*x(3)*(beta1 - 1) + beta1*c*x(3)*(beta2 - 1)) - x(1)*(c*(x(2) - 1) + alpha1*x(2)*(b - c))));
dx(3) = -x(3)*((b*x(3)*(x(1) - 1) - b*beta1*x(1)*x(2))*(x(3) - 1) + (c*(x(3) - 1) + alpha2*x(3)*(b - c))*(x(1) - 1) - x(3)*((c*(x(3) - 1) + alpha2*x(3)*(b - c))*(x(1) - 1) - x(1)*(beta2*c*(x(2) - 1) + beta1*beta2*x(2)*(b - c) - b*beta1*x(2)*(beta2 - 1) + beta2*c*x(2)*(beta1 - 1))) - x(1)*(beta2*c*(x(2) - 1) + beta1*beta2*x(2)*(b - c) - b*beta1*x(2)*(beta2 - 1) + beta2*c*x(2)*(beta1 - 1)));
end

