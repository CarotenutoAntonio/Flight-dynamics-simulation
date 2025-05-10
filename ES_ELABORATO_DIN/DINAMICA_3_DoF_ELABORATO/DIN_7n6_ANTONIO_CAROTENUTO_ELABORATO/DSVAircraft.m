%--------------------------------------------------------------------
% Dinamica e Simulazione di Volo
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
        AR_W
        CL_alpha_W
        eps_0
        DepsDalpha
        %------------------------------------------------------------------
        %Mass, inertia, etc
        %------------------------------------------------------------------
        mass
        W
        k_y
        mu_x
        Xcg_adim
        Xn_adim
        %------------------------------------------------------------------
        %Aerodynamics
        %------------------------------------------------------------------
        CD_0
        K
        m
        CL_alpha
        CL_delta_e
        CL_delta_s
        CL_alpha_dot
        CL_q
        Cm_0
        Cm_alpha
        Cm_alpha_dot
        Cm_delta_s
        Cm_delta_e
        Cm_delta_e_dot
        Ch_e_delta_tab
        Cm_q
        %------------------------------------------------------------------
        %Elevator
        %------------------------------------------------------------------
        S_e
        x_C_e
        Lambda_e
        mac_e
        mass_e
        W_e
        ec_adim
        k_e
        I_e
        I_ey
        Ch_e_0
        Ch_e_alpha
        Ch_e_delta_s
        Ch_e_delta_e
        Ch_e_delta_e_dot
        Ch_e_q
        Ch_e_alpha_dot
        Rs_e
        Rg_e
        Kh_friction
        delta_e_max
        delta_e_min
        %------------------------------------------------------------------
        %Propulsion
        %------------------------------------------------------------------
        T
        Cm_T_0
        Cm_T_alpha
        mu_T
        %------------------------------------------------------------------
        % Limitations
        %------------------------------------------------------------------
        CL_max
        CL_min
        n_max
        n_min
        Fe_max
        Fe_min
    end
    
    methods
        
        %%Constructor, populate class properties reading from file
        function obj = DSVAircraft(dataFileName)
        f_id = fopen(dataFileName,'r');

        %Checking file opening
        if (f_id==-1)
            obj.err = -1; %Opening failed
            disp(['Non è possibile leggere il file.', ...
                dataFileName, ' ...'])
        else
            disp(['File ', ...
                dataFileName, ' aperto.'])
            for i=1:5
                temp = fgetl(f_id); %Read five rows
            end
            
            %%Geometric data
            obj.S = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.b = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.mac = fscanf(f_id,'%f');temp = fgetl(f_id);
            obj.AR_W = ((obj.b)^2)/obj.S;
            for i=1:2
                temp = fgetl(f_id);
            end
            
            %%Mass data
            obj.mass = fscanf(f_id,'%f %*s\n'); temp = fgetl(f_id);
            obj.W = obj.mass*obj.g;
            obj.k_y = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.mu_x = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Xcg_adim = fscanf(f_id,'%f'); temp = fgetl(f_id);
            for i=1:3
                temp = fgetl(f_id);
            end
            
            %%Aerodynamics
            %Neutral point of the aircraft, non-dimensional
            obj.Xn_adim = fscanf(f_id,'%f'); temp = fgetl(f_id);
            temp = fgetl(f_id);
            %Misc
            obj.CL_alpha_W = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.eps_0 = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.DepsDalpha = fscanf(f_id,'%f'); temp = fgetl(f_id);
            temp = fgetl(f_id);
            %Polar
            obj.CD_0 = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.K = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.m = fscanf(f_id,'%f'); temp = fgetl(f_id);
            for i=1:1
                temp = fgetl(f_id);
            end
            %Aerodynamic derivatives
            obj.CL_alpha = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.CL_delta_e = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.CL_delta_s = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.CL_alpha_dot = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.CL_q = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Cm_0 = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Cm_delta_e = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Cm_delta_s = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Cm_alpha_dot = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Cm_q = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Cm_delta_e_dot = fscanf(f_id,'%f'); temp = fgetl(f_id);
            %Stability derivative
            obj.Cm_alpha = -obj.CL_alpha*(obj.Xn_adim - obj.Xcg_adim);
            for i=1:3
                temp = fgetl(f_id);
            end
            
            %%Elevator data
            %Geometry
            obj.S_e = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Lambda_e = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.x_C_e = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.mac_e = fscanf(f_id,'%f');temp = fgetl(f_id);
            temp = fgetl(f_id);
            %Mass, inertia, etc
            obj.mass_e = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.W_e = obj.mass_e*obj.g;
            obj.ec_adim = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.k_e = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.I_e = (obj.mass_e)*(obj.k_e^2);
            obj.I_ey = obj.mass_e*obj.ec_adim*obj.mac_e*obj.x_C_e ...
                - obj.I_e*cos(obj.Lambda_e);
            temp = fgetl(f_id);
            %Aerodynamics
            obj.Ch_e_0 = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Ch_e_alpha = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Ch_e_delta_s = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Ch_e_delta_e = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Ch_e_delta_tab = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Ch_e_delta_e_dot = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Ch_e_q = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Ch_e_alpha_dot = fscanf(f_id,'%f'); temp = fgetl(f_id);
            for i=1:4
                temp = fgetl(f_id);
            end
            
            %%Command linkage characteristics
            %Irreversible Type
            obj.Rs_e = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Rg_e = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Kh_friction = fscanf(f_id,'%f'); temp = fgetl(f_id);
            temp = fgetl(f_id);
            %Angular excursion limitations
            obj.delta_e_max = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.delta_e_min = fscanf(f_id,'%f'); temp = fgetl(f_id);
            for i=1:2
                temp = fgetl(f_id);
            end
            
            %%Propulsion data
            obj.T = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.T = obj.T*obj.g;
            obj.Cm_T_0 = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Cm_T_alpha = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.mu_T = fscanf(f_id,'%f'); temp = fgetl(f_id);
            for i=1:2
                temp = fgetl(f_id);
            end
            
            %%Aerodynamic and structural limitations
            obj.CL_max = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.CL_min = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.n_max = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.n_min = fscanf(f_id,'%f'); temp = fgetl(f_id);
            for i=1:2
                temp = fgetl(f_id);
            end
            
            %%Limitations on piloting stick force
            obj.Fe_max = fscanf(f_id,'%f'); temp = fgetl(f_id);
            obj.Fe_min = fscanf(f_id,'%f'); temp = fgetl(f_id);
            
            %%Setting the error tag
            obj.err = 0;
        end 
        end
        
    end
end