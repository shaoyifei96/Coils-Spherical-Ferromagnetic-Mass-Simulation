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
    end
    methods
        function obj = Field3D(size_map,step_size,mu,Coils,masses)
            if ~mod(size_map,2)%even do
                size_map=size_map+1
            end
            
            grid = -size_map/2:step_size:size_map/2;
            obj.n = length(grid);
            [obj.X,obj.Y,obj.Z] = meshgrid(grid);
            obj.size_map = size_map;
            obj.step_size = step_size;
            obj.mu = mu;
            obj.Coils = Coils;
            obj.masses= masses;
        end
        
        function obj = combineB(obj)
            disp('combineB ran')
            obj.B_total=zeros(obj.n,obj.n,obj.n,3);
            for i=1:length(obj.Coils)
                obj.Coils(i)=obj.Coils(i).Calc_B(obj);
            end
            for i=1:length(obj.Coils)
                obj.B_total(:,:,:,1) = obj.B_total(:,:,:,1)+obj.Coils(i).B(:,:,:,1);%i dir
                obj.B_total(:,:,:,2) = obj.B_total(:,:,:,2)+obj.Coils(i).B(:,:,:,2);%j dir
                obj.B_total(:,:,:,3) = obj.B_total(:,:,:,3)+obj.Coils(i).B(:,:,:,3);%k dir
            end
            Hmagsquared = sum(obj.B_total.^2,4)/(obj.mu^2);%force on mu of enviroment
            [Fx,Fy,Fz] = gradient(Hmagsquared,obj.step_size,obj.step_size,obj.step_size);
            obj.F_total = cat(4,Fx,Fy,Fz);
            
        end
        
        function  [t, state]= simulate_mass(obj,mass_num,t_total,t_precision,movie)
            t=0:t_precision:t_total;
            figure();
            mag2dbcoloredquiver(obj.B_total,obj.Coils,obj,'B total');
            figure();
            mag2dbcoloredquiver(obj.F_total,obj.Coils,obj,'F total');
            xlim([-obj.size_map/2-1 obj.size_map/2+1])
            ylim([-obj.size_map/2-1 obj.size_map/2+1])
            zlim([-obj.size_map/2-1 obj.size_map/2+1])
            mass=obj.masses(mass_num);
            initi_state=[mass.Location mass.Velocity];
            const=2*pi*mass.R^3;%this has problem, missing K, see Jones
            drag_coeff = 4/3*pi*obj.masses(mass_num).R^3*50;%missing Mr %10kg /s/A/m since submerged in fluid may be higher
            [t, state]=ode45(@(t,state)sys(obj,t,state,mass.m,const,drag_coeff),t,initi_state,odeset('RelTol',100,'AbsTol',100));
            
            if movie
                f=figure();
                mag2dbcoloredquiver(obj.F_total,obj.Coils,obj,'F total');
                xlim([-obj.size_map/2-1 obj.size_map/2+1])
                ylim([-obj.size_map/2-1 obj.size_map/2+1])
                zlim([-obj.size_map/2-1 obj.size_map/2+1])
                hold on;
                coeff=floor(length(t)/240);
                h = animatedline('LineWidth',3,'Color','r');
                for j = 1:240
                    addpoints(h,state(j*coeff,1),state(j*coeff,2),state(j*coeff,3));
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
            
            F=obj.Coils.Calc_F_single(obj,reshape(state(1:3),[1 3]),0.01);
            F=const*F;%based on different units may be different, use SI here
            %state=[x  y  z x_dot y_dot z_dot]
            state_dot(1)=state(4);
            state_dot(2)=state(5);
            state_dot(3)=state(6);
            
            force = F(1)-drag_coeff*state(4);
            %             end
            state_dot(4)=force/m;
            %             if abs(drag_coeff*state(5))>1.01*abs(F(2))
            %                 force = -0.01*F(2);
            %             else
            force = F(2)-drag_coeff*state(5);
            %             end
            state_dot(5)=force/m;
            %             if abs(drag_coeff*state(6))>1.01*abs(F(3))
            %                 force = -0.01*F(3);
            %             else
            force = F(3)-drag_coeff*state(6);
            %             end
            state_dot(6)=force/m;
            
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
