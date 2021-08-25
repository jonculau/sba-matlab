classdef LinearSystem < handle

    properties (SetAccess = private)
        y
        x
        A
        B
        C
        D
        Ts
    end
    
    methods
        function obj = LinearSystem(A, B, C, D, Ts, Initials)
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.D = D;
            obj.Ts = Ts;
            obj.x = [];
            obj.y = [];
            if exist('Initials', 'var')  
                obj.x(:, 1) = Initials;
            else
                obj.x(:, 1) = zeros(size(A,1),1);
            end
        end
        
        function y = Sys_Dynamics(obj, u, oldValues, newValues)
            k = size(obj.x, 2);
            
            
            Av = double(subs(obj.A, oldValues, newValues));
            Bv = double(subs(obj.B, oldValues, newValues));
            Cv = double(subs(obj.C, oldValues, newValues));
            Dv = double(subs(obj.D, oldValues, newValues));
            
            if size(Bv,1) ~= 0 
                dx = Av * obj.x(:, k) + Bv * u;
            else
                dx = Av * obj.x(:, k); 
            end
            
            obj.x(:, k + 1) =  dx * obj.Ts + obj.x(:, k);
            
            if size(Dv,1) ~= 0 
                y = Cv * obj.x(:, k + 1) + Dv * u;
                
            else
                y = Cv * obj.x(:, k + 1);
            end
            if k == 1
                obj.y(:,1) = Cv * obj.x(:, 1);
            end
            obj.y(:, k + 1) = y;

        end
    end
end

