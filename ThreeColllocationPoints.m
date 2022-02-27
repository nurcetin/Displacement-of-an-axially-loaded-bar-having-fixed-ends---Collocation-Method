% Collocation Method with Three Collocation Point
clc;
clear;
syms x c1 c2 c3 A ;
% Elasticity Modulus (Pa)
E = 200*10^9;
% Cross Section Area of the Bar (m^2)
A = 3*10^4+1*10^4*x;
% Test Functions
phizero = 0;
phione = x*(2-x);
phitwo = x*x*(2-x);
phithree = x*x*x*(2-x);
% Displacement Function
u = phizero+c1*phione+c2*phitwo+c3*phithree;
% Governing Equation
eqn=diff(E*A*(diff(u)))+3*x;
% Collocation Points
x1=2/4;
x2=2*2/4;
x3=3*2/4;
% Weighted Residual Method
eqn=subs(eqn,x1)==0;
eqn2=subs(eqn,x2)==0;
eqn3=subs(eqn,x3)==0;
% Finding Coefficients c1 and c2
sol = solve([eqn1, eqn2, eqn3], [c1, c2, c3]);
c1 = sol.c1;
c2 = sol.c2;
c3 = sol.c3;
% Approximate Solution for Displacement u(x)
u = phizero+c1*phione+c2*phitwo+c3*phithree;
% Approximate Solution for N(x)
N = E*A*diff(u);
% Check for Boundary Conditions
boundary_conditions=subs(u,0);
boundary_conditions=subs(u,2);
