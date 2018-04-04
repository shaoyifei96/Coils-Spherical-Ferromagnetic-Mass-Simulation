function f = EllipticE(z,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: D. Cebron
% MATLAB implementation of the following Maple function
% EllipticE - Incomplete and complete elliptic integrals of the second kind
% z = algebraic expression (the sine of the amplitude)
% k = algebraic expression (the parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(length(varargin)==1),
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Incomplete elliptic integrals of the second kind
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k=varargin{1};
    
    % MATLAB package Elliptic_Integrals.zip by Thomas Hoffend
    f=lellipe(asin(z),k,1e-13);
    
    % MATLAB package Elliptic Integrals and Functions by Moiseev Igor
    % [~,f,~] = elliptic12(asin(z),k^2); % if asin(z) is real
    % or
    % [f,~,~] = elliptic12i(asin(z),k^2) % if asin(z) is complex
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Complete elliptic integrals of the second kind
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % MATLAB package Elliptic_Integrals.zip by Thomas Hoffend
    f=lellipe(pi/2,z,1e-13);
    
    % MATLAB package Elliptic Integrals and Functions by Moiseev Igor
    % [~,f,~] = elliptic12(pi/2,z^2);
end