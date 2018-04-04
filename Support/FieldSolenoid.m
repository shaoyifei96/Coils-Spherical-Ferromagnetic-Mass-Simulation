function B=FieldSolenoid(a,L,Br,rho,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: D. Cebron
% Field at (rho,z) in polar coordinates of a cylindrical magnet
% (modeled as a solenoid) radius a, length L, remanent flux density Br
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some useful quantities for the calculation
zp=z+L./2; zm=z-L./2;
kp=sqrt(4.*a.*rho./((a+rho).^2+zp.^2)); km=sqrt(4.*a.*rho./((a+rho).^2+zm.^2));

%%%%%%%%%%%%%%%%%%% RADIAL FIELD OF THE SOLENOID %%%%%%%%%%%%%%%%%%%%%%%%%%
% Callaghan 1960 with correction of the sign error
% (or Wikipedia with k in argument of Elliptic functions)

if(rho==0) % Field on the solenoid axis
    B(1)=0;
else
    B(1)=Br./pi.*sqrt(a./rho).*(((kp.^2-2)./(2.*kp).*EllipticK(kp)+...
        EllipticE(kp)./kp)-((km.^2-2)./(2.*km).*EllipticK(km)+EllipticE(km)./km));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% AXIAL FIELD OF THE SOLENOID %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derby Olbert 2009
% (after a MAPLE conversion in terms of Elliptic integrals)
Bo=Br./pi;
kp=sqrt((zp.^2+(a-rho).^2)./(zp.^2+(a+rho).^2));
km=sqrt((zm.^2+(a-rho).^2)./(zm.^2+(a+rho).^2));
betap=zp./sqrt(zp.^2+(rho+a).^2); betam=zm./sqrt(zm.^2+(rho+a).^2);
g=(a-rho)./(a+rho);

if(kp<1)
    B(2)=Bo.*a./(a+rho).*1./(g+1).*(betap.*(EllipticK(sqrt(1-kp.^2))+g.*EllipticPi(1-g.^2,sqrt(1-kp.^2)))...
        -betam.*(EllipticK(sqrt(1-km.^2))+g.*EllipticPi(1-g.^2,sqrt(1-km.^2))));
else
    B(2)=Bo.*a./(a+rho).*1./(g.*(g+1)).*(betap./kp.*(EllipticK(sqrt((-1+kp).*(kp+1))./kp).*g+...
        EllipticPi((-1+g).*(g+1)./g.^2, sqrt((-1+kp).*(kp+1))./kp))...
        -betam./km.*(EllipticK(sqrt((-1+km).*(km+1))./km).*g+EllipticPi((-1+g).*(g+1)./g.^2, sqrt((-1+km).*(km+1))./km)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%