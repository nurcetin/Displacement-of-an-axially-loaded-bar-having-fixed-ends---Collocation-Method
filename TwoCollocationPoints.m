% Collocation Method with Two Collocation Points
clc;
clear;
syms x c1 c2 A ;
% Elasticity Modulus (Pa)
E = 200*10^9;
% Cross Section Area of the Bar (m^2)
A = 3*10^4+1*10^4*x;
% Test Functions
phizero = 0;
phione = x*(2-x);
phitwo = x*x*(2-x);
% Displacement Function
u = phizero+c1*phione+c2*phitwo;
% Governing Equation
eqn=diff(E*A*(diff(u)))+3*x;
% Collocation Points
x1=2/3;
x2=2*2/3;
% Weighted Residual Method
eqn=subs(eqn,x1)==0;
eqn2=subs(eqn,x2)==0;
% Finding Coefficients c1 and c2
sol = solve([eqn1, eqn2], [c1, c2]);
c1 = sol.c1;
c2 = sol.c2;
% Approximate Solution for Displacement u(x)
u = phizero+c1*phione+c2*phitwo;
% Approximate Solution for N(x)
N = E*A*diff(u);
% Check for Boundary Conditions
boundary_conditions=subs(u,0);
boundary_conditions=subs(u,2);
