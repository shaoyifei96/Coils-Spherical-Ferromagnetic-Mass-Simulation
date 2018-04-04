% Make the hgtransform and the surface; parent the surface to the hgtransform
[x,y,z] = cylinder(1);
h = hgtransform;
surf(x,y,z, 'Parent', h, 'FaceColor', 'r');
view(3)
% Make it taller
% Tip it over and make it taller
set(h, 'Matrix', makehgtform('xrotate', 1, 'scale', [1 1 1],'translate',1,1,1))
% ax = axes('XLim',[-1.5 1.5],'YLim',[-1.5 1.5],'ZLim',[-1.5 1.5]);
% view(3)
% grid on
% 
% [x,y,z] = cylinder([.2 0]);
% h(1) = surface(x,y,z,'FaceColor','red');
% h(2) = surface(x,y,-z,'FaceColor','green');
% h(3) = surface(z,x,y,'FaceColor','blue');
% h(4) = surface(-z,x,y,'FaceColor','cyan');
% h(5) = surface(y,z,x,'FaceColor','magenta');
% h(6) = surface(y,-z,x,'FaceColor','yellow');
% 
% t = hgtransform('Parent',ax);
% set(h,'Parent',t)
% 
% Rz = eye(4);
% Sxy = Rz;
% 
% for r = 1:.1:2*pi
%     % Z-axis rotation matrix
%     Rz = makehgtform('zrotate',r);
%     % Scaling matrix
%     Sxy = makehgtform('scale',r/4);
%     % Concatenate the transforms and
%     % set the transform Matrix property
%     set(t,'Matrix',Rz*Sxy)
%     drawnow
% end