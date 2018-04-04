function f = EllipticF(z,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: D. Cebron
% MATLAB implementation of the following Maple function
% EllipticF - Incomplete elliptic integral of the first kind
% z = algebraic expression (the sine of the amplitude)
% k = algebraic expression (the parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MATLAB package Elliptic_Integrals.zip by Thomas Hoffend
f=lellipf(asin(z),k,1e-13);

% MATLAB package Elliptic Integrals and Functions by Moiseev Igor
%[f,~,~] = elliptic12(asin(z),k^2); % if asin(z) is real
% or
%[f,~,~] = elliptic12i(asin(z),k^2) % if asin(z) is complex
