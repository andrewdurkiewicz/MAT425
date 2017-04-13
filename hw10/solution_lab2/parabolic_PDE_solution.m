%% Solve the elliptic PDE
%%    ∂_t u = α^2 ∂_x^2 u    ,  0<x<l, t>0
%%       u(x,0) = u0(x)      ,  0<x<l, t=0
%%     u(0,t) = u(l,t) = 0   ,     t≥0


%%-- parameters
alpha2 = 2;
l      = 10;
T      = 1;
u0 = @(x) (x>=4).*(x<=6);

%%-- numerical paramters
dx = .5;
dt = .005;
ld = alpha2/dx^2*dt;

%% init
nTime = floor(T/dt+.5);
nX    = floor(l/dx+.5);

%% IC
intX  = linspace(0,l,nX);
u     = u0(intX);


%%----------------------------------%%
%%-------      loop         --------%%
%%----------------------------------%%

## for n=1:nTime
##   %% new time step
##   u_new = u;
##   for i=2:(nX-1)
##     u_new(i) = (1-2*ld)*u(i) + ld*(u(i-1)+u(i+1));
##   end
##   %% plot
##   plot(intX,u_new)
##   xlabel('x')
##   axis([0 l 0 1])
##   title(['solution at time t=' num2str(n*dt,'%1.2f')])
##   %% update
##   u = u_new;
##   pause(.0001)
## end


for n=1:nTime
  %% new time step
  u = (1-2*ld)*u + ld*([u(2:end) 0] + [0 u(1:(end-1))]);
  trapz(intX,u)
  %% plot
  plot(intX,u)
  xlabel('x')
  axis([0 l 0 1])
  title(['solution at time t=' num2str(n*dt,'%1.2f')])
  %% update
  pause(.0001)
end



