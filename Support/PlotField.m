%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: D. Cebron
% Show how to plot magnetic fields using the scripts 
% FieldSolenoid & FieldBar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1 solenoid
clear all; close all; clc;
a=1;        % Solenoid radius
L=1;        % Solenoid length
Br=1;       % Residual field = field in infinite solenoid (mu*n*I)

[X,Y] = meshgrid(-4*max(a,L):max(a,L)/(3*pi):4*max(a,L));
for k=1:length(X)
    for kk=1:length(Y)
        rho=X(k,kk); z=Y(k,kk);
        BB=FieldSolenoid(a,L,Br,rho,z); % BB(1:2)=(Brho,Bz)
        Brho=BB(1); Bz(k,kk)=BB(2);
        Bx(k,kk)=Brho*sign(X(k,kk)); 
        B(k,kk)=sqrt(BB(1)^2+BB(2)^2);
    end
end
figure(1); contour(X,Y,B,100,'EdgeColor','none'); colorbar; hold on; %colormap('Hot')
plot([-1 1]*a,[1 1]*L/2,'-k','Linewidth',2); plot([-1 1]*a,[1 1]*(-L/2),'-k','Linewidth',2);
plot([1 1]*a,[-1 1]*L/2,'-k','Linewidth',2); plot(-[1 1]*a,[-1 1]*(L/2),'-k','Linewidth',2);
quiver(X(1:5:end,1:5:end),Y(1:5:end,1:5:end),Brho(1:5:end,1:5:end),Bz(1:5:end,1:5:end),4,'g','Linewidth',2);
h = streamline(X,Y,Bx,Bz,X(1:10:end,1:10:end),Y(1:10:end,1:10:end)); set(h,'Color','red');
%% 2 solenoid
clear all; clc;
a=1;        % Solenoid radius
L=1;        % Solenoid length
Br=1;       % Residual field  = field in infinite solenoid (mu*n*I)
Wz=2*a;     % Intermagnet distance
hz=Wz+L;    % Distance between magnets centers

[X,Y] = meshgrid(-4*max(a,L):max(a,L)/(3*pi):4*max(a,L));
for k=1:length(X)
    for kk=1:length(Y)
        rho=X(k,kk); z=Y(k,kk)+(Wz+L)/2;
        BB=FieldSolenoid(a,L,Br,rho,z)...
            +FieldSolenoid(a,L,Br,rho,z-hz); % BB(1:2)=(Brho,Bz)
        Brho=BB(1); Bz(k,kk)=BB(2);
        Bx(k,kk)=Brho*sign(X(k,kk)); Bz(k,kk)=BB(2); 
        B(k,kk)=sqrt(BB(1)^2+BB(2)^2);
    end
end
figure(2); contourf(X,Y,B,100,'EdgeColor','none'); colorbar; hold on; %colormap('Hot')
plot([-1 1]*a,[1 1]*L/2-(Wz+L)/2,'-k','Linewidth',2); plot([-1 1]*a,[1 1]*(-L/2)-(Wz+L)/2,'-k','Linewidth',2);
plot([1 1]*a,[-1 1]*L/2-(Wz+L)/2,'-k','Linewidth',2); plot(-[1 1]*a,[-1 1]*(L/2)-(Wz+L)/2,'-k','Linewidth',2);
plot([-1 1]*a,[1 1]*L/2+hz-(Wz+L)/2,'-k','Linewidth',2); plot([-1 1]*a,[1 1]*(-L/2)+hz-(Wz+L)/2,'-k','Linewidth',2);
plot([1 1]*a,[-1 1]*L/2+hz-(Wz+L)/2,'-k','Linewidth',2); plot(-[1 1]*a,[-1 1]*(L/2)+hz-(Wz+L)/2,'-k','Linewidth',2);
% quiver(X,Y,Brho,Bz,4,'g');
h = streamline(X,Y,Bx,Bz,X(1:25:end,1:25:end),Y(1:25:end,1:25:end)); set(h,'Color','red');
%% 1 cuboidal magnet
clear all; clc;
Lz=1;       % Magnet size in the z-direction
Lx=16*Lz;   % Magnet size in the x-direction 
Ly=3.4*Lz;  % Magnet size in the y-direction
Br=1;       % Residual field 

[X,Z] = meshgrid(-15*Lz:Lz/(3*pi):15*Lz);
for k=1:length(X)
    for kk=1:length(Z)
        BB=FieldBar(Lx,Ly,Lz,Br,X(k,kk),0,Z(k,kk));
        Bx(k,kk)=BB(1); By(k,kk)=BB(2); Bz(k,kk)=BB(3); B(k,kk)=sqrt(BB(1)^2+BB(3)^2);
    end
end
figure(3); contourf(X,Z,B,100,'EdgeColor','none'); colorbar; hold on; %colormap('Hot')
plot([-1 1]*Lx/2,[1 1]*Lz/2,'-k','Linewidth',2); plot([-1 1]*Lx/2,[1 1]*(-Lz/2),'-k','Linewidth',2);
plot([1 1]*Lx/2,[-1 1]*Lz/2,'-k','Linewidth',2); plot(-[1 1]*Lx/2,[-1 1]*(Lz/2),'-k','Linewidth',2);
% quiver(X(1:3:end,1:3:end),Z(1:25:end,1:25:end),Bx(1:25:end,1:25:end),Bz(1:25:end,1:25:end),1,'k');
h = streamline(X,Z,Bx,Bz,X(1:25:end,1:25:end),Z(1:25:end,1:25:end)); set(h,'Color','red');
%% 2 cuboidal attractive magnets
clear all; clc;
Lz=1;       % Magnet size in the z-direction
Lx=16*Lz;   % Magnet size in the x-direction 
Ly=3.4*Lz;  % Magnet size in the y-direction
Br=1;       % Residual field 
Wz=6*Lz;    % Intermagnet distance
hz=Wz+Lz;   % Distance between magnets centers

[X,Z] = meshgrid(-15*Lz:Lz/(3*pi):15*Lz);
for k=1:length(X)
    for kk=1:length(Z)
        BB=FieldBar(Lx,Ly,Lz,Br,X(k,kk),0,Z(k,kk)+(Wz+Lz)/2)+FieldBar(Lx,Ly,Lz,Br,X(k,kk),0,Z(k,kk)-hz+(Wz+Lz)/2);
        Bx(k,kk)=BB(1); By(k,kk)=BB(2); Bz(k,kk)=BB(3); B(k,kk)=sqrt(BB(1)^2+BB(3)^2); 
    end
end
figure(4);
%contourf(X,Z,B,100,'EdgeColor','none'); colorbar; hold on; %colormap('Hot')
plot([-1 1]*Lx/2,[1 1]*Lz/2-(Wz+Lz)/2,'-k','Linewidth',2); plot([-1 1]*Lx/2,[1 1]*(-Lz/2)-(Wz+Lz)/2,'-k','Linewidth',2);
plot([1 1]*Lx/2,[-1 1]*Lz/2-(Wz+Lz)/2,'-k','Linewidth',2); plot(-[1 1]*Lx/2,[-1 1]*(Lz/2)-(Wz+Lz)/2,'-k','Linewidth',2);
plot([-1 1]*Lx/2,[1 1]*Lz/2+hz-(Wz+Lz)/2,'-k','Linewidth',2); plot([-1 1]*Lx/2,[1 1]*(-Lz/2)+hz-(Wz+Lz)/2,'-k','Linewidth',2);
plot([1 1]*Lx/2,[-1 1]*Lz/2+hz-(Wz+Lz)/2,'-k','Linewidth',2); plot(-[1 1]*Lx/2,[-1 1]*(Lz/2)+hz-(Wz+Lz)/2,'-k','Linewidth',2);
% quiver(X(1:3:end,1:3:end),Z(1:3:end,1:3:end),Bx(1:3:end,1:3:end),Bz(1:3:end,1:3:end),4,'k');
h = streamline(X,Z,Bx,Bz,X(1:25:end,1:25:end),Z(1:25:end,1:25:end)); set(h,'Color','red');