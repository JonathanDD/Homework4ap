close all;clear all;clc;

%% Discretization
% Create amount to discretize by

N=2^8*10;

% discretize and plug in variables

A=1;

k=100;
v=1;
x = linspace(0,1,N);
a = -2*k^2*(1/N)^2*ones(1,N);
b = ones(1,N);
c = ones(1,N);
c(1)=2;
u = zeros(1,N);
f = A*ones(1,N);
alpha = zeros(1,N);

%% Apply the boundary conditions


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
    uxact(l) = [1-cosh(k*x(l))/cosh(k*x(N))]*A/k^2-v/k*sinh(k*(x(N)-x(l)))/cosh(k*x(l));
end

uofy = @(y) [1-[sinh(10*(1-y))+sinh(10*y)]/sinh(10*y)]*1/10^2+1*sinh(10*(1-y))/sinh(10*y);

subplot(2,1,1)
plot(x,u)
subplot(2,1,2)
plot(x,uxact)
%subplot(3,1,3)
%fplot(uofy, [ 0 1])