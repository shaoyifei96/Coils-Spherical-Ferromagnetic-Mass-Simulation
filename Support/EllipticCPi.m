function f = EllipticCPi(nu,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: D. Cebron
% MATLAB implementation of the following Maple function
% EllipticCPi - Complementary complete elliptic integral of the third kind
% nu =algebraic expression (the characteristic)
% k = algebraic expression (the parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = EllipticPi(nu,sqrt(1-k.^2));