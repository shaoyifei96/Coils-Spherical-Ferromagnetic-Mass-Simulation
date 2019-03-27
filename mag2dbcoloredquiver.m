function mag2dbcoloredquiver(B_total,Coils,field,tile,twod)
ax=axes;
n=field.n;
Bmag   =zeros(n,n,n);
B_unit =zeros(n,n,n,3);
if twod
    for i =1:n
        for j =1:n
            k=(n-1)/2+1
                Bmag(i,j,k) = sqrt(sum(B_total(i,j,k,:).^2));
                B_unit(i,j,k,:)=B_total(i,j,k,:)./Bmag(i,j,k);

                %B_unit_Radial(i,j,k,:)=reshape(Radial_unit(i,j,k,:),[1 3])/sqrt(sum(reshape(Radial_unit(i,j,k,:),[1 3]).^2));
        end
    end
else
    for i =1:n
        for j =1:n
            for k =1:n
                Bmag(i,j,k) = sqrt(sum(B_total(i,j,k,:).^2));
                B_unit(i,j,k,:)=B_total(i,j,k,:)./Bmag(i,j,k);
                %B_unit_Radial(i,j,k,:)=reshape(Radial_unit(i,j,k,:),[1 3])/sqrt(sum(reshape(Radial_unit(i,j,k,:),[1 3]).^2));
                %B_unit_Normal(i,j,k,:)=obj.Direction;
            end
        end
    end
end
q=quiver3(field.X,field.Y,field.Z,B_unit(:,:,:,1),B_unit(:,:,:,2),B_unit(:,:,:,3));
%// Create a quiver3 as we normally would (could also be 2D quiver)


%// Compute the magnitude of the vectors
%mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
%    reshape(q.WData, numel(q.UData), [])).^2, 2));
Bmag(isnan(Bmag))=1e-10;
Bmag=(mag2db(Bmag));
mags = reshape(Bmag,[n.^3 1]);
%// Get the current colormap
currentColormap = colormap(hsv);

%// Now determine the color to make each arrow using a colormap
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));

%// Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'

%// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');
%quiver3(field.X,field.Y,field.Z,B_unit(:,:,:,1),B_unit(:,:,:,2),B_unit(:,:,:,3));
%quiverC3D(field.X(10,:,:),field.Y(10,:,:),field.Z(10,:,:),obj.B(10,:,:,1),obj.B(10,:,:,2),obj.B(10,:,:,3),500);
%slice(field.X,field.Y,field.Z,mag2db(obj.Bmag),[2.5],[2.5],[0:5])
xlabel('x');
ylabel('y');
zlabel('z');
colorbar;
x0=100;
y0=100;
width=1000;
height=800;
set(gcf,'units','points','position',[x0,y0,width,height]);
title(tile);

cylinders    = zeros(length(Coils));
target_tf    = zeros(length(Coils));
hold on;
for i= 1:length(Coils)

[x_cy,y_cy,z_cy]= cylinder(Coils(i).R,40);
z_cy=z_cy-0.5;
cylinders(i)=surf(x_cy,y_cy,z_cy, 'FaceColor', 'r');
target_tf(i)   = hgtransform('Parent',ax);
set(cylinders(i),'Parent', target_tf(i));
[theta,phi,~] = cart2sph(Coils(i).Direction(1),Coils(i).Direction(2),Coils(i).Direction(3));
translate_form = makehgtform('translate',Coils(i).Location(1),Coils(i).Location(2),Coils(i).Location(3));%+Coils(i).L/2);
scale_form = makehgtform('scale', [1 1 Coils(i).L]);
rotate_form = makehgtform('zrotate',theta,'yrotate',-pi/2+phi);
set(target_tf(i),'Matrix',translate_form*rotate_form*scale_form );%
%
end
hold off;
axis equal;

end

