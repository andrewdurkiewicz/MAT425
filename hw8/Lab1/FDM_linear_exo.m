%% Solve the BVP with the finite difference method
%%   y'' = p(x) y' + q(x) y + r(x)
%%   y(a) = α, y(b) = β

%% Boundary value problem
a     = 0;
b     = 1;
alpha = 2;
beta  = 2;
p = @(x) -1./(x+1);
q = @(x) 2 + 0*x;
r = @(x) (1-x.^2).*exp(-x);

%% Parameters
N  = 10-1;
dx = (b-a)/(N+1);
x  = a+dx*(1:N);
% explanation
%------------
% 0      1      2                          N     N+1
% |------|------|------| ... |------|------|------|
% a     a+dx   a+2dx                              b


%% Initialization
A = ...
vecB = -dx^2*r(x');
%% boundary condition
vecB(1) = ...
vecB(N) = ...


%%-----       solve        ---%%
%%----------------------------%%
Y = A\vecB;
% add the boundary condition
y_sol = [alpha; Y; beta];


%%----- plot
plot([a x b],y_sol')
xlabel('x')
ylabel('y')
title('y''''=p y + q y'' + r, y(a)=alpha,y(b)=beta')
