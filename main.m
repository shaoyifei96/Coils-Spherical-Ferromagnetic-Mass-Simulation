clear all; clc; close all;
addpath('Support');

I=3; %Current A
N=5000;% nuber of turns
L=0.08;%length m
R=0.01;% Radius m

mu=1.2566370614e-6;% enviroment permeability
map_size= 0.1;%quiverdb plot has to have odd sized map size.
%if a even number is entered, it automatically add 1 to it
%if even, also print the actual size in the command line
twod = 1;
%if want 2d, toggle twod to 1
step_size= 0.023;%map resolution;for visualization purpose only; don't use a number that can be divided by map size

Location1=[-0.06 0 0];
Direction1=[1 0 0];

Location2=[0 0.06 0];
Direction2=[0 -1 0];

Location3=[1 5 0];
Direction3=[0 -1 0];

%Create Coils by calling the calss 
%ElectroMagnet_Class(0:dont use eatimation, Current, Turns, Length, Radius, Location Vector3, Direction Vector3)
Coil1 = ElectroMagnet_Class(0,I  ,N,L,R,Location1,Direction1);
Coil2 = ElectroMagnet_Class(0,I/6,N,L,R,Location2,Direction2);
Coil3 = ElectroMagnet_Class(0,I,N,0.1,R*10,Location3,Direction3);

%combine into a structure
Coils(1)=Coil1;
Coils(2)=Coil2;
%Coils(3)=Coil3;


mu2=1e3;%1e3 -1e5
sus=50; %20-70

R_seph=0.001;
m=0.0001;%kg
Location4=[0 0 0];
Velocity4=[0 0 0];

%ferro_linear_mass(Radius,Permeability,mass,Location vector3, Direction Vector3);
mass1=ferro_linear_mass(R_seph,mu2,m,Location4, Velocity4,sus);

R_seph=0.01;
m=0.002;%kg
Location4=[0 0 0];
Velocity4=[0 0 0];

%mass2=ferro_linear_mass(R_seph,mu2,m,Location4, Velocity4);
 
%Create 3D field with certain size and [array of Coils] and ONE mass
%Field3D(map_size,step_size,permeability,[array of Coils],mass object,2d_switch);
Field1 = Field3D(map_size,step_size,mu,Coils,mass1,twod);

Field1=Field1.combineB();%this is important
%this step calculate the linearly superimporsed B field
%if don't do this, there is no B field anywhere, you would encounter
%an error on mag2dbcoloredquiver of index exceeds dimension.

%Field1.simulate_mass(total time length(s),each time step length(s),movie?);
[t, state]= Field1.simulate_mass(300,0.06,1);%when closing a animation when not done, output a error with invalid/deleted obj.
%when using simulation, if animation does not work corrently, try changing
%the time step(try both increase and decrease time step)
% % % % hold on;
% % % %             scatter3(state(:,1),state(:,2),state(:,3));
% % % %             p = scatter3(state(1,1),state(1,2),state(1,3),150,'green','filled');
% % % %             hold off;
% % % %             for k = 2:length(t)
% % % %                 p.XData = state(k,1)
% % % %                 p.YData = state(k,2);
% % % %                 p.ZData = state(k,3);
% % % %                 drawnow
% % % %             end
% % % %         end

%other useful functions

%to visualize a vector field(with magnitude represented by color)
%use mag2dbcoloredquiver(vector field,[array of coils](can be empty),Field_obj(required),Title);

%to observe some peoperty follow the following structure:
%Field3D:
% 	size_map
% 	step_size
% 	mu
% 	Coils[ElectroMagnet_Class ]
% 		ElectroMagnet_Class
% 			Diopole_est bool
% 			I
% 			N
% 			L
% 			R
% 			Location
% 			Direction
% 			B(from just 1 Coil)
% 	masses[ferro_linear_mass] %currently support only 1 mass
% 		ferro_linear_mass
% 			R
% 			mu
% 			m
% 			Location
% 			Direction
% 	X meshgrid
% 	Y meshgrid
% 	Z meshgrid
% 	B_total : 4D vector
% 				dim1: x
% 				dim2: y
% 				dim3: z 
% 				dim4:i j k
% 	F_total : 4D vector
% 				dim1: x
% 				dim2: y
% 				dim3: z 
% 				dim4:i j k

% To get current of Coil1: enter in command line Field1.Coils(1).I





