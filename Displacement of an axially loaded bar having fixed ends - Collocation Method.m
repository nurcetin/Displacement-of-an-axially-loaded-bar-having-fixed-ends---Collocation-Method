% Collocation Method with One Collocation Point
clc;
clear;
syms x c1 A ;
% Elasticity Modulus (Pa)
E = 200*10^9;
% Cross Section Area of the Bar (m^2)
A = 3*10^4+1*10^4*x;
% Test Functions
phizero = 0;
phione = x*(2-x);
% Displacement Function
u = phizero+c1*phione;
% Governing Equation
eqn=diff(E*A*(diff(u)))+3*x;
% Collocation Point x1
x1=2/2;
% Weighted Residual Method
eqn=subs(eqn,x1)==0;
% Finding Coefficient c1
c1 = solve(eqn,c1);
% Approximate Solution for Displacement u(x)
u = phizero+c1*phione;
% Approximate Solution for N(x)
N = E*A*diff(u);
% Check for Boundary Conditions
boundary_conditions=subs(u,0);
boundary_conditions=subs(u,2);

