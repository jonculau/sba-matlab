classdef StewartPlatform < handle
    properties        
        % System Parameters
        
        m = 1.36                                          % Body mass
        Im = diag([1.705e-4;1.705e-4; 3.404e-4])          % Tensor of inertia
        rB = 0.35;                                        % Bottom radius;
        phiB = pi / 6;                                    % Angle of each leg to bottom;
        eaB = Quaternion.ea2q([0 0 0]);                                   % Bottom orientation;
        rT = 0.25;                                        % Top radius;
        phiT = pi / 2;                                    % Angle of each leg to top;        
        
        
        % Home Position
        tB = [0; % 'x' home position [m];
            0; % 'y' home position [m];
            0]; % 'z' home position [m];
        
        tT = [0; % 'x' home position [m];
            0; % 'y' home position [m];
            0.57145]; % 'z' home position [m];
    end
    properties (SetAccess = private)
        % States regarding the top platform local frame
        x = [];
        
        %Inputs
        u = [];           % System inputs
        d = [];
        
        dot_l = [];       % Linear actuators force
        dot_ld = [];
        %Simulation Parameters
        Ts % time stamp
        
        B;      % Bottom platform vectors
        
        T;      % Top platform vectors
        
        L;      % Leg vectorss
        
    end
    properties (Dependent)
        J = zeros(6, 12); % Jacobian Matrix
        eaT = [];           % Euler angles conversion
        t                   % Time vector
    end
    
    methods
%% Get Methods
        function [J] = get.J(obj)
            J = zeros(6, 6);
            k = size(obj.x, 2);
            q = obj.x(k).q;
            p = obj.x(k).p;
            for i = 1:6
                L = (q.toR * obj.T(1:3, i) + p) - obj.B(1:3, i);
                n = L./sqrt(sum(L.^2));
                J(i, :) = [n' (skew(q.toR * obj.T(1:3, i)) * n)']; 
            end
        end        
        
        function [eaT] = get.eaT(obj)             
            eaT = zeros(3,size(obj.x,2));
            for i = 1:size(obj.x,2)
                eaT(:,i) = obj.x(i).q.toEulerAngle;
            end
        end
        
        function [t] = get.t(obj)
            t = 0:obj.Ts:size(obj.x,2)*obj.Ts - obj.Ts;
        end
        %% Simulation Methods
       
        % Reset Parameters
        function obj = StewartPlatform(Ts, InitialCondtions)
            obj.x = [];
            if(~exist('InitialCondtions', 'var'))
                obj.x(1).q = Quaternion([1;0;0;0]);
                obj.x(1).w = zeros(3,1);
                obj.x(1).p = obj.tT;
                obj.x(1).v = zeros(3,1);
            else
                if class(InitialCondtions.q) ~= "Quaternion" 
                    error("The initial value of q must be a quaternion")
                end
                obj.x(1).q = InitialCondtions.q;
                obj.x(1).w = InitialCondtions.w;
                obj.x(1).p = InitialCondtions.p;
                obj.x(1).v = InitialCondtions.v;
            end
            if (~exist('ts', 'var'))
                 obj.Ts = 1e-3;
            else
                 obj.Ts = Ts;
            end
            
            obj.u = [];
            obj.dot_l = [];
            
            [obj.B,obj.T, obj.L] = obj.Make;
        end
        
        % Sys dynamics
        function [x] = Sys_Dynamics(obj, inputs, disturbance)
            
            k = size(obj.x, 2); %Current indice
            
            u = obj.J' * (inputs + disturbance);
            
            % Storage system input values
            obj.dot_l(:, k) = inputs;  
            obj.dot_ld(:, k) = disturbance; 
            obj.u (:,k) =  obj.J' * (inputs);
            obj.d (:,k) =  obj.J' * (disturbance);
            
            obj.dot_l(:, k + 1) = inputs; % Para plot gráfico            
            obj.u (:,k + 1) =  obj.J' * inputs; % Para plot gráfico 
            obj.d (:,k + 1) =  obj.J' * (disturbance); % Para plot gráfico 
            obj.dot_ld(:, k + 1) = disturbance; % Para plot gráfico 
            
            % Get current values
            x = obj.x(k);
            
            % Computate derivatives            
            dq = zeros(4, 1);
            
            dq(1) = x.q.e' * x.w;
            
            dq(2:4) = (x.q.n * eye(3) + skew(x.q.e)) * x.w/2;
            
            dw = obj.Im \ (u(1:3) - skew(x.w) * obj.Im * x.w); 
            
            dv = obj.m \ u(4:6);
            
            dp = x.v;

            % Storage new values
            obj.x(k + 1).q = Quaternion([x.q.n;x.q.e] + obj.Ts * dq);
            
            obj.x(k + 1).w = dw * obj.Ts + x.w;
            
            obj.x(k + 1).v = dv * obj.Ts + x.v;
            
            obj.x(k + 1).p = dp * obj.Ts + x.p;
            
            x = obj.x(k + 1);
            

        end


        function AnimationPlot(obj)
            px = 0; % 'x' home position [m];
            py = 0; % 'y' home position [m];
            pz = 0.57145 - (1.7/1000); % 'z' home position [m];
            
            % Environment:
            figure();
            axis3plot([0 0 0]', [0 0 0]', '-', 1.5);
            view([150 20]);
            offset = .8;
            axis([px - offset px + offset py - offset py + offset 0 - 0.2 pz + 0.2]);
            title('Pose Estimation');
            ylabel('y'); xlabel('x'); zlabel('z');
            grid on;
            
            % Stewart:
            
            % Poses
            gT = [obj.x(1).q.toR	obj.x(1).p; 0	0	0	1] * obj.T;
            
            gB = [obj.eaB.toR obj.tB; 0 0 0 1] * obj.B;
            
            plotB = plot3([gB(1, :) gB(1, 1)]', [gB(2, :) gB(2, 1)]', [gB(3, :) gB(3, 1)]', 'k', 'Linewidth', 2);
            plotT = patch(gT(1, :)', gT(2, :)', gT(3, :)', 'y', 'FaceAlpha', 0);
            plotS(6) = 0;
            
            for j = 1:6
                S = [gB(:, j) gT(:, j)];
                plotS(j) = plot3(S(1, :)', S(2, :)', S(3, :)', 'r--', 'Linewidth', 1);
            end
            
            % Animation:
            disp('Animating...')
            
            for i = 1:10:size(obj.x,2)
                % Poses:
                gT = [
                    obj.x(i).q.toR	obj.x(i).p
                    0	0	0	1
                    ];
                
                gT = gT * obj.T;
                set(plotT, 'XData', gT(1, :), 'YData', gT(2, :), 'ZData', gT(3, :));
                
                % Stems:
                for j = 1:6
                    S = [gB(:, j) gT(:, j)];
                    set(plotS(j), 'XData', S(1, :), 'YData', S(2, :), 'ZData', S(3, :));
                end
                
                drawnow;
            end
            
        end


    end

  methods (Access = private)      
      % Graphical functions
      function [B, T, L] = Make(this)
          
          % Joints Definitions (Local Frame):
          lambda = zeros(1, 6); % Angle of each vertice;
          v = [];
          
          for i = 1:2:5
              lambda(i) = i * pi / 3 - this.phiB / 2;
              lambda(i + 1) = lambda(i) + this.phiB;
              
              v(i) = i * pi / 3 - this.phiT / 2;
              v(i + 1) = v(i) + this.phiT;
          end
          
          B = zeros(4, 6); % Point of each vertice (bottom);
          T = zeros(4, 6); % Point of each vertice (top);
          L = zeros(1, 6);
          
          for i = 1:6
              B(:, i) = [this.rB * cos(lambda(i)) this.rB * sin(lambda(i)) 0 1]';
              T(:, i) = [this.rT * cos(v(i)) this.rT * sin(v(i)) 0 1]';
              L(i) = sqrt((B(:, i) - T(:, i))' * (B(:, i) - T(:, i)));
          end
          
      end
  end

end
