% This class defines Electromagnets. Requires  packages: eliptical Intergal to compute
% the Biot Savat. Doesn't need any other package if using Dipole Estimate.
%                 Dipole_est: 0
%             I: 5
%             N: 2000
%             L: 0.2000
%             R: 0.6000
%      Location: [-5 0 0]
%     Direction: [0.5774 0.5774 0.5774] %auto convert to unit vector
%             B: []
% By Yifei Shao 2018
classdef ElectroMagnet_Class
    properties (Access = public)
        Dipole_est
        I
        N
        L
        R
        Location
        Direction
        
        B%the 4dimension one with i j k
    end
    
    properties (SetAccess = private, GetAccess = ? Field3D)
        Br
        Bn
        Bmag
        
    end
    methods
        %constructor
        function obj = ElectroMagnet_Class(Dipole_est,I,N,L,R,Location,Direction)
            obj.I = I;
            obj.N = N;
            obj.L = L;
            obj.R = R;
            obj.Dipole_est=Dipole_est;
            if all(size(Location) ==[1 3]) && all(size(Direction) ==[1 3])
                obj.Location  = Location;
                obj.Direction = Direction./norm(Direction);
            else
                error('please use [x y z] for Coil Location and Direction');
            end
        end
        % Determine the Difference bewteen a BiotSavat Calc and Point
        % Dipole Estimation and Graph it
        
        function Graph_Error(obj,field)
            Brr=obj.N*obj.I*field.mu;  %mu0     % Residual field = field in infinite solenoid (mu*n*I)
            n=length(field.X);
            B_mag=zeros(n,n);
            for i=1:n
                for j=1:n
                    BB=FieldSolenoid(obj.R,obj.L,Brr,field.X(i,j,1),field.Y(i,j,1)); % BB(1:2)=(Brho,Bz)
                    B_mag(i,j)=sqrt(BB(1)^2+BB(2)^2);
                end
            end% compute accurate intergral
            
            Bx = field.mu*(obj.I*(obj.R^2)/4).*((field.X(:,:,1).^2+field.Y(:,:,1).^2).^(-3/2)).*(2-(3.*field.Y(:,:,1).^2)./(field.X(:,:,1).^2+field.Y(:,:,1).^2));%mu0
            By = field.mu*((3*obj.I*(obj.R^2)/4).*field.X(:,:,1).*field.Y(:,:,1))./((field.X(:,:,1).^2+field.Y(:,:,1).^2).^(5/2));%mu0
            Bmag_est=sqrt(Bx.^2+By.^2);
            Est_diff=abs(Bmag_est-B_mag);%Compute Dipole Estimate
            
            surf(field.X(:,:,1),field.Y(:,:,1),mag2db(Est_diff));
            colormap winter
            xlabel('x');
            ylabel('y');
            zlabel('dB|B|')
            leg=strcat('Coil L=',num2str(obj.L),' on y axis');
            legend(leg)
            x0=100;
            y0=100;
            width=1000;
            height=800;
            set(gcf,'units','points','position',[x0,y0,width,height])%Formatting
            
        end
        
        %calculate Br Bn and B magnitude
        function F=Calc_F_single(obj,field,loc,precision)% this can be for multiple coils together
            %get the location of the point with respect to the Coil
            x_low =loc-[precision 0 0];
            x_high=loc+[precision 0 0];
            y_low =loc-[0 precision 0];
            y_high=loc+[0 precision 0];
            z_low =loc-[0 0 precision];
            z_high=loc+[0 0 precision];
            loc=cat(1,x_low,x_high,y_low,y_high,z_low,z_high);
            B_point = zeros(length(loc),3);
            for n=1:length(obj)
                for i=1:length(loc)
                    Dist_Vec=loc(i,:)-obj(n).Location;
                    Normal_Dist=dot(Dist_Vec,obj(n).Direction);
                    Normal_Vec=Normal_Dist * obj(n).Direction;
                    
                    Radial_Vec= Dist_Vec-Normal_Vec;
                    Radial_Dist=sqrt(sum(Radial_Vec.^2));
                    Radial_unit=Radial_Vec./Radial_Dist;
                    if isnan(Radial_unit)
                        Radial_unit=[1 0 0];
                    end
                    Brr=obj(n).N*obj(n).I*field.mu;
                    
                    BB=FieldSolenoid(obj(n).R,obj(n).L,Brr,Radial_Dist,Normal_Dist); % BB(1:2)=(Brho,Bz)
                    Br_single=BB(1); Bn_single=BB(2);
                    B_point(i,:)  = B_point(i,:)+Radial_unit.* Br_single+ obj(n).Direction.*Bn_single;
                end
            end
            % calculate H^2
            B_point=sum(B_point.^2,2);
            Hmagsq=B_point/(field.mu^2);
            d=precision*2;
            F=[(Hmagsq(2)-Hmagsq(1))/d (Hmagsq(4)-Hmagsq(3))/d (Hmagsq(6)-Hmagsq(5))/d];
            
            
        end
        
        function obj=Calc_B(obj,field)
            %Let z axis be the up
            n=length(field.X);
            Normal_Vec=zeros(n,n,n,3);
            Normal_Dist=zeros(n,n,n);
            Dist_Vec=zeros(n,n,n,3);
            
            %             X_vec = field.X-obj.Location(1);
            %             Y_vec = field.Y-obj.Location(2);
            %             Z_vec = field.Z-obj.Location(3);
            for i =1:n
                for j =1:n
                    for k =1:n
                        %optimize the following by not using for loop
                        Normal_Dist(i,j,k)=dot(...
                            [field.X(i,j,k)-obj.Location(1)...
                            field.Y(i,j,k)-obj.Location(2)...
                            field.Z(i,j,k)-obj.Location(3)],obj.Direction);
                        Normal_Vec(i,j,k,:)=Normal_Dist(i,j,k) * obj.Direction;
                        
                        Dist_Vec(i,j,k,:)=[field.X(i,j,k) field.Y(i,j,k) field.Z(i,j,k)] - obj.Location;
                        
                    end
                end
            end
            Radial_Vec= Dist_Vec-Normal_Vec;
            Radial_Dist=sqrt(sum(Radial_Vec.^2, 4));
            Radial_unit=Radial_Vec./Radial_Dist;
            
            %Normal Vector not necessary
            if obj.Dipole_est %use dipole estimation
                r=(Radial_Dist.^2+Normal_Dist.^2).^0.5;
                obj.Br = (field.mu*obj.I*(obj.R^2)/4).*((r).^(-3)).*(2-(3.*Normal_Dist.^2)./(r.^2));
                obj.Bn = ((3*field.mu*obj.I*(obj.R^2)/4).*Radial_Dist.*Normal_Dist)./(r.^5);
                
            else%when estimation is off, use biot savat
                
                Brr=obj.N*obj.I*field.mu;
                for i=1:n
                    for j=1:n
                        for k=1:n
                            BB=FieldSolenoid(obj.R,obj.L,Brr,Radial_Dist(i,j,k),Normal_Dist(i,j,k)); % BB(1:2)=(Brho,Bz)
                            obj.Br(i,j,k)=BB(1); obj.Bn(i,j,k)=BB(2);
                        end
                    end
                end% compute accurate intergral
            end
            
            
                          obj.B=zeros(n,n,n,3);
            %             B_unit=zeros(n,n,n,3);
            for i =1:n
                for j =1:n
                    for k =1:n
                        if any(isnan(Radial_unit(i,j,k,:)))
                        Radial_unit(i,j,k,:)=[1 0 0];
                        end
                        obj.B(i,j,k,:)  = reshape(Radial_unit(i,j,k,:),[1 3]).* obj.Br(i,j,k)+ obj.Direction.*obj.Bn(i,j,k);
                        %                         obj.Bmag(i,j,k) = sqrt(sum(obj.B(i,j,k,:).^2));
                        %                         B_unit(i,j,k,:)=obj.B(i,j,k,:)./obj.Bmag(i,j,k);
                        %B_unit_Radial(i,j,k,:)=reshape(Radial_unit(i,j,k,:),[1 3])/sqrt(sum(reshape(Radial_unit(i,j,k,:),[1 3]).^2));
                        %B_unit_Normal(i,j,k,:)=obj.Direction;
                    end
                end
            end

        end
        
        
         function obj=Calc_B_2d(obj,field)
            %Let z axis be the up
            n=length(field.X);
            Normal_Vec=zeros(n,n,n,3);
            Normal_Dist=zeros(n,n,n);
            Dist_Vec=zeros(n,n,n,3);
            
            %             X_vec = field.X-obj.Location(1);
            %             Y_vec = field.Y-obj.Location(2);
            %             Z_vec = field.Z-obj.Location(3);
            for i =1:n
                for j =1:n
                    k= (n-1)/2+1;
                        %optimize the following by not using for loop
                        Normal_Dist(i,j,k)=dot(...
                            [field.X(i,j,k)-obj.Location(1)...
                            field.Y(i,j,k)-obj.Location(2)...
                            0],obj.Direction);
                        Normal_Vec(i,j,k,:)=Normal_Dist(i,j,k) * obj.Direction;
                        
                        Dist_Vec(i,j,k,:)=[field.X(i,j,k) field.Y(i,j,k) field.Z(i,j,k)] - obj.Location;
                        
                    
                end
            end
            Radial_Vec= Dist_Vec-Normal_Vec;
            Radial_Dist=sqrt(sum(Radial_Vec.^2, 4));
            Radial_unit=Radial_Vec./Radial_Dist;
            
             
             
                Brr=obj.N*obj.I*field.mu;
                for i=1:n
                    for j=1:n
                        k= (n-1)/2+1;
                            BB=FieldSolenoid(obj.R,obj.L,Brr,Radial_Dist(i,j,k),Normal_Dist(i,j,k)); % BB(1:2)=(Brho,Bz)
                            obj.Br(i,j,k)=BB(1); obj.Bn(i,j,k)=BB(2);
                    end
                end% compute accurate intergra
            
                         obj.B=zeros(n,n,n,3);
            %             B_unit=zeros(n,n,n,3);
            for i =1:n
                for j =1:n
                    k= (n-1)/2+1;
                        if any(isnan(Radial_unit(i,j,k,:)))
                        Radial_unit(i,j,k,:)=[1 0 0];
                        end
                        obj.B(i,j,k,:)  = reshape(Radial_unit(i,j,k,:),[1 3]).* obj.Br(i,j,k)+ obj.Direction.*obj.Bn(i,j,k);
                        %                         obj.Bmag(i,j,k) = sqrt(sum(obj.B(i,j,k,:).^2));
                        %                         B_unit(i,j,k,:)=obj.B(i,j,k,:)./obj.Bmag(i,j,k);
                        %B_unit_Radial(i,j,k,:)=reshape(Radial_unit(i,j,k,:),[1 3])/sqrt(sum(reshape(Radial_unit(i,j,k,:),[1 3]).^2));
                        %B_unit_Normal(i,j,k,:)=obj.Direction;
                end
            end

        end
        
    end
end

            %mag2dbcoloredquiver(obj.B,obj,field,strcat('B ',' Location=',num2str(obj.Location),' Direction=',num2str(obj.Direction)));
            %             q=quiver3(field.X,field.Y,field.Z,B_unit(:,:,:,1),B_unit(:,:,:,2),B_unit(:,:,:,3));
            %             %// Create a quiver3 as we normally would (could also be 2D quiver)
            %
            %
            %             %// Compute the magnitude of the vectors
            %             %mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
            %             %    reshape(q.WData, numel(q.UData), [])).^2, 2));
            %             obj.Bmag(isnan(obj.Bmag))=0;
            %             obj.Bmag=(mag2db(obj.Bmag));
            %             mags = reshape(obj.Bmag,[field.n.^3 1]);
            %             %// Get the current colormap
            %             currentColormap = colormap(hot);
            %
            %             %// Now determine the color to make each arrow using a colormap
            %             [~, ~, ind] = histcounts(mags, size(currentColormap, 1));
            %
            %             %// Now map this to a colormap to get RGB
            %             cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
            %             cmap(:,:,4) = 255;
            %             cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
            %
            %             %// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
            %             set(q.Head, ...
            %                 'ColorBinding', 'interpolated', ...
            %                 'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
            %
            %             %// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
            %             set(q.Tail, ...
            %                 'ColorBinding', 'interpolated', ...
            %                 'ColorData', reshape(cmap(1:2,:,:), [], 4).');
            %             %quiver3(field.X,field.Y,field.Z,B_unit(:,:,:,1),B_unit(:,:,:,2),B_unit(:,:,:,3));
            %             %quiverC3D(field.X(10,:,:),field.Y(10,:,:),field.Z(10,:,:),obj.B(10,:,:,1),obj.B(10,:,:,2),obj.B(10,:,:,3),500);
            %             %slice(field.X,field.Y,field.Z,mag2db(obj.Bmag),[2.5],[2.5],[0:5])
            %             xlabel('x');
            %             ylabel('y');
            %             zlabel('z');
            %             colorbar;
            %             leg1=strcat('Coil Location =',num2str(obj.Location),'Coil Direction = ',num2str(obj.Direction));
            %             legend(leg1,'Location','best');
            %
            %             return;