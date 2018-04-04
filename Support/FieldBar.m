function B=FieldBar(Lx,Ly,Lz,Br,x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: D. Cebron
% Magnetic flux density B at (x,y,z) created by a bar magnet 
% of size (Lx,Ly,Lz) and uniform magnetization M0=mu*Br along the z axis
% (coordinate centered on the magnet center, oriented by its inertia axes)
% From Camacho & Sosa Rev. Mex. E, E 59, 8-17, 2013
% (B in T if input data are S.I.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu0=4*pi*1e-7;
M0=Br/mu0;
a=Lx/2; b=Ly/2; c=Lz/2; % (Lx,Ly,Lz)=(2*a,2*b,2*c)

% From Camacho & Sosa Rev. Mex. E, E 59, 8-17, 2013
% (after a MAPLE calculation of the expressions)
B(1)=1/4*mu0*M0/pi*log(((x^2-2*x*a+a^2+y^2-2*y*b+b^2+z^2-2*z*c+c^2)^(1/2)+b-y)/((x^2-2*x*a+a^2+y^2+2*y*b+b^2+z^2-2*z*c+c^2)^(1/2)-b-y)*((x^2+2*x*a+a^2+y^2-2*y*b+b^2+z^2+2*z*c+c^2)^(1/2)+b-y)/((x^2+2*x*a+a^2+y^2+2*y*b+b^2+z^2+2*z*c+c^2)^(1/2)-b-y)/((x^2+2*x*a+a^2+y^2-2*y*b+b^2+z^2-2*z*c+c^2)^(1/2)+b-y)*((x^2+2*x*a+a^2+y^2+2*y*b+b^2+z^2-2*z*c+c^2)^(1/2)-b-y)/((x^2-2*x*a+a^2+y^2-2*y*b+b^2+z^2+2*z*c+c^2)^(1/2)+b-y)*((x^2-2*x*a+a^2+y^2+2*y*b+b^2+z^2+2*z*c+c^2)^(1/2)-b-y));
B(2)=1/4*mu0*M0/pi*log(((y^2-2*y*a+a^2+x^2-2*x*b+b^2+z^2-2*z*c+c^2)^(1/2)+b-x)/((y^2-2*y*a+a^2+x^2+2*x*b+b^2+z^2-2*z*c+c^2)^(1/2)-b-x)*((y^2+2*y*a+a^2+x^2-2*x*b+b^2+z^2+2*z*c+c^2)^(1/2)+b-x)/((y^2+2*y*a+a^2+x^2+2*x*b+b^2+z^2+2*z*c+c^2)^(1/2)-b-x)/((y^2+2*y*a+a^2+x^2-2*x*b+b^2+z^2-2*z*c+c^2)^(1/2)+b-x)*((y^2+2*y*a+a^2+x^2+2*x*b+b^2+z^2-2*z*c+c^2)^(1/2)-b-x)/((y^2-2*y*a+a^2+x^2-2*x*b+b^2+z^2+2*z*c+c^2)^(1/2)+b-x)*((y^2-2*y*a+a^2+x^2+2*x*b+b^2+z^2+2*z*c+c^2)^(1/2)-b-x));
Bz=-1/4*mu0*M0/pi*(atan((-x+a)*(y+b)/(z+c)/(x^2-2*x*a+a^2+y^2+2*y*b+b^2+z^2+2*z*c+c^2)^(1/2))+atan((-x+a)*(y+b)/(-z+c)/(x^2-2*x*a+a^2+y^2+2*y*b+b^2+z^2-2*z*c+c^2)^(1/2))+atan((-x+a)*(-y+b)/(z+c)/(x^2-2*x*a+a^2+y^2-2*y*b+b^2+z^2+2*z*c+c^2)^(1/2))+atan((-x+a)*(-y+b)/(-z+c)/(x^2-2*x*a+a^2+y^2-2*y*b+b^2+z^2-2*z*c+c^2)^(1/2))+atan((x+a)*(y+b)/(z+c)/(x^2+2*x*a+a^2+y^2+2*y*b+b^2+z^2+2*z*c+c^2)^(1/2))+atan((x+a)*(y+b)/(-z+c)/(x^2+2*x*a+a^2+y^2+2*y*b+b^2+z^2-2*z*c+c^2)^(1/2))+atan((x+a)*(-y+b)/(z+c)/(x^2+2*x*a+a^2+y^2-2*y*b+b^2+z^2+2*z*c+c^2)^(1/2))+atan((x+a)*(-y+b)/(-z+c)/(x^2+2*x*a+a^2+y^2-2*y*b+b^2+z^2-2*z*c+c^2)^(1/2)));

if((abs(x)<=a)&&(abs(y)<=b)&&(abs(z)<=c))
    B(3)=Bz+mu0*M0;
else
    B(3)=Bz;
end