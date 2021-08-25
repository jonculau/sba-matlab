
classdef Quaternion 
    properties (SetAccess = private)
        n
        e
    end
    methods
        function obj = Quaternion(V)
            n = sqrt(sum(V.^2));
            q = V./ (ones([size(V,1) 1])* n);
            obj.n = q(1);
            obj.e = q(2:4);
        end  
        function [ea] = toEulerAngle(obj)
            ea = quat2eul([obj.n obj.e'])';
        end
        function R = toR(obj)
            q1 = obj.n;
            q2 = obj.e(1);
            q3 = obj.e(2);
            q4 = obj.e(3);            
            R = [ q1^2+q2^2-q3^2-q4^2  2*(q2*q3-q1*q4)      2*(q2*q4+q1*q3)
            2*(q2*q3+q1*q4)      q1^2-q2^2+q3^2-q4^2  2*(q3*q4-q1*q2)
            2*(q2*q4-q1*q3)      2*(q3*q4+q1*q2)      q1^2-q2^2-q3^2+q4^2];
        end
    end
    methods (Static)
        function [q] = ea2q(ea)
            if size(ea,2) == 3
                ea = ea';
            end
            cang = cos( ea/2 );
            sang = sin( ea/2 );

            q = Quaternion([-sang(1,:).*sang(2,:).*sang(3,:) + cang(1,:).*cang(2,:).*cang(3,:) 
            sang(1,:).*cang(2,:).*cang(3,:) + sang(2,:).*sang(3,:).*cang(1,:)
            -sang(1,:).*sang(3,:).*cang(2,:) + sang(2,:).*cang(1,:).*cang(3,:) 
            + sang(1,:).*sang(2,:).*cang(3,:) + sang(3,:).*cang(1,:).*cang(2,:)]);
        
        end 
    end
end

