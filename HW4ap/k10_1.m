close all;clear all;clc;

%% Discretization
% Create amount to discretize by

N=2^5*10;

% discretize

A=1;

k=10;
x = linspace(0,1,N);
a = -2*k^2*(1/N)^2*ones(1,N);
b = ones(1,N);
c = ones(1,N);
u = zeros(1,N);
f = A*ones(1,N);
alpha = zeros(1,N);

%% Apply the boundary conditions

u(1)= 10;
u(N)=0;

alpha(1)=a(1);
g(1)=f(1);

for j=2:1:N
    alpha(j)=a(j)-(b(j)/alpha(j-1))*c(j-1);
    g(j)=f(j) - (b(j)/alpha(j-1))*g(j-1);
end


for k = 1:1:(N-1)
    u(N-k)= (g(N-k)-c(N-k)*u(N-k))/alpha(N-k);
end

%% Solve uxact for comparison

for l=1:N
    uxact(l) = [1-[sinh(k*(x(N)-x(l)))+sinh(k*x(l))]/sinh(k*x(l))]*A/k^2+10*sinh(k*(x(N)-x(l)))/sinh(k*x(l));
end

subplot(2,1,1)
plot(x,u)
subplot(2,1,2)
plot(x,uxact)