classdef Field3D
    properties (Access = public)
        size_map
        step_size
        mu
        n
        Coils
        masses
        
        X
        Y
        Z
        
        B_total
        F_total
        Hmagsquared
        
        twod
    end
    methods
        function obj = Field3D(size_map,step_size,mu,Coils,masses,twod)
            
            grid = 0:step_size:(size_map+0.1)/2;
            grid = [fliplr(-grid) grid(2:end)];
            %grid = -size_map/2:step_size:size_map/2;
            obj.n = length(grid);
            [obj.X,obj.Y,obj.Z] = meshgrid(grid);
            obj.size_map = size_map;
            obj.step_size = step_size;
            obj.mu = mu;
            obj.Coils = Coils;
            obj.masses= masses;
            obj.twod = twod;
            
        end
        
        function obj = combineB(obj)
            disp('combineB ran')
            obj.B_total=zeros(obj.n,obj.n,obj.n,3);
            if obj.twod
                for i=1:length(obj.Coils)
                    obj.Coils(i)=obj.Coils(i).Calc_B_2d(obj);
                end

            else
            for i=1:length(obj.Coils)
                obj.Coils(i)=obj.Coils(i).Calc_B(obj);
            end     
            end
            for i=1:length(obj.Coils)
                obj.B_total(:,:,:,1) = obj.B_total(:,:,:,1)+obj.Coils(i).B(:,:,:,1);%i dir
                obj.B_total(:,:,:,2) = obj.B_total(:,:,:,2)+obj.Coils(i).B(:,:,:,2);%j dir
                obj.B_total(:,:,:,3) = obj.B_total(:,:,:,3)+obj.Coils(i).B(:,:,:,3);%k dir
            end
            obj.Hmagsquared = sum(obj.B_total.^2,4)/(obj.mu^2);%force on mu of enviroment
            
            [Fx,Fy,Fz] = gradient(obj.Hmagsquared,obj.step_size,obj.step_size,obj.step_size);
            obj.F_total = cat(4,Fx,Fy,Fz);
            
        end
        
        function  [t, state]= simulate_mass(obj,t_total,t_precision,movie)
            t=0:t_precision:t_total;
            figure(1);
            mag2dbcoloredquiver(obj.B_total,obj.Coils,obj,'B total (Color Scale in dB)',obj.twod);
            for i=1 : length(obj.Coils)
                loc=obj.Coils(i).Location;
                text(loc(1),loc(2),loc(3),['N*I = ' num2str(obj.Coils(i).N*obj.Coils(i).I)],'Color','Magenta')
            end
            xlabel('x coordinate (m)');
            ylabel('y coordinate (m)');
            zlabel('z coordinate (m)');
            
            figure(2);
            mag2dbcoloredquiver(obj.F_total,obj.Coils,obj,'F total (Color Scale in dB)',obj.twod);
            for i=1 : length(obj.Coils)
                loc=obj.Coils(i).Location;
                text(loc(1),loc(2),loc(3),['N*I = ' num2str(obj.Coils(i).N*obj.Coils(i).I)],'Color','Magenta');
            end
            xlabel('x coordinate (m)');
            ylabel('y coordinate (m)');
            zlabel('z coordinate (m)');
            xlim([-obj.size_map/2-obj.size_map/10 obj.size_map/2+obj.size_map/10])
            ylim([-obj.size_map/2-obj.size_map/10 obj.size_map/2+obj.size_map/10])
            zlim([-obj.size_map/2-obj.size_map/10 obj.size_map/2+obj.size_map/10])
            state = zeros(length(obj.masses),length(t),6);
            for mass_num = 1:length(obj.masses) 
                mass=obj.masses(mass_num);
                initi_state=[mass.Location mass.Velocity];
                %const=2*pi*mass.R^3/3*obj.mu*mass.sus/(1+mass.sus/3);
                %since eq achieved in nanoseconds, ignore drag and always use ss v (Shipiro) drag_coeff = 4/3*pi*obj.masses(mass_num).R^3*50;%missing Mr %10kg /s/A/m since submerged in fluid may be higher
                [t_output, state(mass_num,:,:)]=ode45(@(t,state)sys(obj,t,state),t,initi_state);odeset('RelTol',1000,'AbsTol',1000)
                
            end
            
            
            if movie
                f=figure(3);
                mag2dbcoloredquiver(obj.F_total,obj.Coils,obj,'Motion of Particle (Quiver:F total (Color Scale in dB))',obj.twod);
                xlim([-obj.size_map/2-obj.size_map/10 obj.size_map/2+obj.size_map/10]);
                ylim([-obj.size_map/2-obj.size_map/10 obj.size_map/2+obj.size_map/10]);
                zlim([-obj.size_map/2-obj.size_map/10 obj.size_map/2+obj.size_map/10]);
                xlabel('x coordinate (m)');
                ylabel('y coordinate (m)');
                zlabel('z coordinate (m)');
                for i=1 : length(obj.Coils)
                    loc=obj.Coils(i).Location;
                    text(loc(1),loc(2),loc(3),['N*I = ' num2str(obj.Coils(i).N*obj.Coils(i).I)],'Color','Magenta');
                end
                
                for mass_num = 1:length(obj.masses) 
                    loc = obj.masses(mass_num).Location;
                    vel = obj.masses(mass_num).Velocity;
                    text(loc(1),loc(2),loc(3),['Loc = ' num2str(loc) newline 'Vel = ' num2str(vel) newline 'Mass = ' num2str(obj.masses(mass_num).m) ],'Color','Magenta');
                    vel_text(mass_num) = text (loc(1),loc(2),loc(3),'-');
                    hold on;
                    coeff=length(t_output)/240;
                    h(mass_num) = animatedline('LineWidth',3,'Color','r');
                end
                if coeff<1
                    frames = length(t_output);
                else
                    frames = 240;
                end
                for j = 1:frames
                    for mass_num = 1:length(obj.masses)     
                        vector_num = floor(j*coeff);
                        if vector_num == 0
                            vector_num = 1;
                        end
                        %vector_num
                        addpoints(h(mass_num),state(mass_num,vector_num,1),state(mass_num,vector_num,2),state(mass_num,vector_num,3));
                        vel_text(mass_num).String = ['v_x = ' num2str(state(mass_num,vector_num,4)) newline 'v_y = ' num2str(state(mass_num,vector_num,5)) newline 'v_z = ' num2str(state(mass_num,vector_num,6))];
                        vel_text(mass_num).Position =[state(mass_num,vector_num,1) state(mass_num,vector_num,2) state(mass_num,vector_num,3)]; 
                    end
                    drawnow
                    Frames(j) = getframe(f);
                    percent=j/240*100
                end
                v = VideoWriter('output.avi');
                v.FrameRate = 24;  % Default 30
                %v.Quality = 80;    % Default 75
                disp("writing to file...");
                open(v);
                writeVideo(v, Frames);
                close(v);
                disp("done");
            end
            
            
            
        end
        
        
        
        function state_dot=sys(obj,t,state,m,const,drag_coeff)
            F=obj.Coils.Calc_F_single(obj,reshape(state(1:3),[1 3]),0.0001);%0.5mm length scale
            v=4e-13*F;%3.5e-13  -  4.2e-13 m^4/A^2 (Probst)
            %state=[x  y  z x_dot y_dot z_dot]
            state_dot(1)=v(1);
            state_dot(2)=v(2);
            state_dot(3)=v(3);
            
            
            state_dot(4)=0;
            %             if abs(drag_coeff*state(5))>1.01*abs(F(2))
            %                 force = -0.01*F(2);
            %             else
            
            state_dot(5)=0;
            %             if abs(drag_coeff*state(6))>1.01*abs(F(3))
            %                 force = -0.01*F(3);
            %             else
%             force = F(3)-drag_coeff*state(6);
%             %             end
            state_dot(6)=0;
            
            state_dot=state_dot(:);
            
            
        end
        
        
        %
        %             obj.Location()
        %
        %             out(2)=1/obj.m*obj.F();
        %
        %             out(4)=1/obj.m*obj.F();
        %
        %             out(6)=1/obj.m*obj.F();
        %         end
        
    end
    
end
