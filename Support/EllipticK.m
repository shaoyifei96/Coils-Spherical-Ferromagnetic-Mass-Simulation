function f = EllipticK(k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: D. Cebron
% MATLAB implementation of the following Maple function
% EllipticK - Complete elliptic integral of the first kind
% z = algebraic expression (the sine of the amplitude)
% k = algebraic expression (the parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MATLAB package Elliptic_Integrals.zip by Thomas Hoffend
f=lellipf(pi/2,k,1e-13);

% MATLAB package Elliptic Integrals and Functions by Moiseev Igor
%[f,~,~] = elliptic12(pi/2,k^2);
