function f = EllipticPi(z,nu,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: D. Cebron
% MATLAB implementation of the following Maple function
% EllipticPi - Incomplete and complete elliptic integrals of the third kind
% z = algebraic expression (the sine of the amplitude)
% nu =algebraic expression (the characteristic)
% k = algebraic expression (the parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(length(varargin)==1),
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Incomplete elliptic integrals of the third kind
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    k=varargin{1};
    % MATLAB package Elliptic_Integrals.zip by Thomas Hoffend
    f=lellippi(asin(z),-k,-nu,1e-13);
    
    % MATLAB package Elliptic Integrals and Functions by Moiseev Igor
    % f=elliptic3(asin(z),k^2,nu);
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Complete elliptic integrals of the third kind
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % MATLAB package Elliptic_Integrals.zip by Thomas Hoffend
    f=lellippi(pi/2,nu,-z,1e-13);
    
    % MATLAB package Elliptic Integrals and Functions by Moiseev Igor
    % f=elliptic3(pi/2,nu^2,z);
end
