function f = EllipticCK(k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: D. Cebron
% MATLAB implementation of the following Maple function
% EllipticCK - Complementary complete elliptic integral of the first kind
% z = algebraic expression (the sine of the amplitude)
% k = algebraic expression (the parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = EllipticK(sqrt(1-k.^2));