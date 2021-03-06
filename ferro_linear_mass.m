% This is just a struct to hold properties for ferromagnetic  mass
% assume no saturation. 
% theory based on Jones Electromechanics of Particles
% Yifei Shao 2018
classdef ferro_linear_mass 
    properties (Access = public)
        R
        mu
        m  
        Location
        Velocity
        sus
    end
    methods
        function obj = ferro_linear_mass(R,mu,m,Location, Velocity,sus)
            obj.R=R;
            obj.m=m;
            obj.mu = mu;
            obj.Location=Location;
            obj.Velocity=Velocity;
            obj.sus = sus; %suseptibility
            %K=(mu-Field.mu)/(mu+2*Field.mu);
        end
        
       
    end
    
end
