%--------------------------------------------------------------------
% Dinamica e Simulazione di Volo - Raia Andrea
%--------------------------------------------------------------------
classdef DSVAircraft
%DSVAircraft class for aircraft data management
%Detailed explanation goes here
   properties
        %------------------------------------------------------------------
        %IDs and misc.
        %------------------------------------------------------------------
        Name = 'DSVAircraft - <Put a name here>';
        g = 9.81;
        err = 0;
        %------------------------------------------------------------------
        %Geometry
        %------------------------------------------------------------------
        S
        b
        mac
        %------------------------------------------------------------------
        %Mass, inertia, etc
        %------------------------------------------------------------------
        mass
        W
        Ixx
        Iyy
        Izz
        Ixy
        Ixz
        Iyz
        I1
        I2
        I3
        I4
        I5
        I6
        I_matrix
        det_I
        %------------------------------------------------------------------
        %Propulsion
        %------------------------------------------------------------------
        T_max_SL
        e_T
        mu_T
    end
    
    methods
        
        
         %% Constructor, populate class properties reading from file
         function obj = DSVAircraft(dataFileName)
         f_id = fopen(dataFileName,'r');
            
         %Checking file opening
         if (f_id==-1)
             obj.err = -1; %Opening failed
             disp(['DSVAircraft :: initFromFile __ Could NOT open file ', ...
                 dataFileName, ' ...'])
         else
             disp(['DSVAircraft :: initFromFile __ Opening file ', ...
                 dataFileName, ' ... OK.'])
             for i=1:5
                 temp = fgetl(f_id); %Read five rows
             end
             
             %%Geometric data
             obj.S = fscanf(f_id,'%f'); 
             obj.S = convlength(convlength(obj.S,'ft','m'),'ft','m'); temp = fgetl(f_id);
             obj.b = fscanf(f_id,'%f');
             obj.b = convlength( obj.b,'ft','m'); temp = fgetl(f_id);
             obj.mac = fscanf(f_id,'%f');
             obj.mac = convlength (obj.mac,'ft','m'); temp = fgetl(f_id);
             for i=1:2
                 temp = fgetl(f_id);
             end
             
             %%Mass data
             obj.mass = fscanf(f_id,'%f %*s\n'); 
             obj.mass = convmass(obj.mass,'slug','kg'); temp = fgetl(f_id);
             obj.W = obj.mass*obj.g;
             obj.Ixx = fscanf(f_id,'%f'); 
             obj.Ixx = convmass(obj.Ixx,'slug','kg');
             obj.Ixx = convlength(convlength(obj.Ixx,'ft','m'),'ft','m'); temp = fgetl(f_id);
             obj.Iyy = fscanf(f_id,'%f'); 
             obj.Iyy = convmass(obj.Iyy,'slug','kg');
             obj.Iyy = convlength(convlength(obj.Iyy,'ft','m'),'ft','m'); temp = fgetl(f_id);     
             obj.Izz = fscanf(f_id,'%f'); 
             obj.Izz = convmass(obj.Izz,'slug','kg');
             obj.Izz = convlength(convlength(obj.Izz,'ft','m'),'ft','m'); temp = fgetl(f_id);              
             obj.Ixy = fscanf(f_id,'%f'); 
             obj.Ixy = convmass(obj.Ixy,'slug','kg');
             obj.Ixy = convlength(convlength(obj.Ixy,'ft','m'),'ft','m'); temp = fgetl(f_id);              
             obj.Ixz = fscanf(f_id,'%f'); 
             obj.Ixz = convmass(obj.Ixz,'slug','kg');
             obj.Ixz = convlength(convlength(obj.Ixz,'ft','m'),'ft','m'); temp = fgetl(f_id);
             obj.Iyz = fscanf(f_id,'%f'); 
             obj.Iyz = convmass(obj.Iyz,'slug','kg');
             obj.Iyz = convlength(convlength(obj.Iyz,'ft','m'),'ft','m'); temp = fgetl(f_id);
             obj.I1 = obj.Iyy*obj.Izz - (obj.Iyz)^2;
             obj.I2 = obj.Ixy*obj.Izz + obj.Iyz*obj.Ixz;
             obj.I3 = obj.Ixy*obj.Iyz + obj.Iyy*obj.Ixz;
             obj.I4 = obj.Ixx*obj.Izz - (obj.Ixz)^2;
             obj.I5 = obj.Ixx*obj.Iyz + obj.Ixy*obj.Ixz;
             obj.I6 = obj.Ixx*obj.Iyy - (obj.Ixy)^2;
             obj.I_matrix = [obj.Ixx   -obj.Ixy   -obj.Ixz;
                             -obj.Ixy   obj.Iyy   -obj.Iyz;
                             -obj.Ixz  -obj.Iyz    obj.Izz];
             obj.det_I = (obj.Ixx*obj.Iyy*obj.Izz)-...
                         (obj.Ixx*(obj.Iyz)^2)-...
                         (obj.Izz*(obj.Ixy)^2)-...
                         (obj.Iyy*(obj.Ixz)^2)-...
                         (2*obj.Iyz*obj.Ixz*obj.Ixy);
             for i=1:2
                 temp = fgetl(f_id);
             end
             
             %%Propulsion data
             obj.T_max_SL = fscanf(f_id,'%f');
             obj.T_max_SL = convforce(obj.T_max_SL,'lbf','N'); temp = fgetl(f_id);
             obj.e_T = fscanf(f_id,'%f'); 
             obj.e_T = convlength(obj.e_T,'ft','m'); temp = fgetl(f_id);
             obj.mu_T = fscanf(f_id,'%f'); temp = fgetl(f_id);
                
             %%Setting the error tag
             obj.err = 0;
         end
         end
    end
end