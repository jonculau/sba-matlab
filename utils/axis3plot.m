function  axis3plot(o,R,linha,width)
%axis3plot=[o,D] Plots an axis gimble at coordinate o=[x y z] 
%with direction D = [teta,fi,psi]

d = 50; %length of axis

s = min(size(R));

if s== 1
    t = R(1);
    f = R(2);
    p = R(3);

  %rotation matrices
  Rt = [1  0       0
        0  cos(t) -sin(t)
        0  sin(t)  cos(t)];
  Rf = [cos(f) 0 sin(f)
        0      1 0
        -sin(f) 0 cos(f)];     
  Rp = [cos(p) -sin(p) 0 
        sin(p) cos(p) 0
        0         0    1];
  R  = Rt*Rf*Rp;
end


%%
if nargin == 2    
    linha = '-';
    width = 1;
elseif nargin ==3
    width = 1;    
elseif nargin < 2
    disp('Not enought arguments.')
    disp('Input must be o=[x,y,z] and D = [teta,fi,psi]')
    return
elseif nargin > 5
    disp('Too many arguments.')
    disp('Input must be o=[x,y,z] and D = [teta,fi,psi]')
    return
end

%%
o = [o(1) o(2) o(3)]';


%compute axes
X= R*[d 0 0]' + o;
Y= R*[0 d 0]' + o;
Z= R*[0 0 d]' + o;

plot3([o(1) X(1)],[o(2) X(2)],[o(3) X(3)],strcat('g',linha),'Linewidth',width) %x
hold on
plot3([o(1) Y(1)],[o(2) Y(2)],[o(3) Y(3)],strcat('r',linha),'Linewidth',width) %y
plot3([o(1) Z(1)],[o(2) Z(2)],[o(3) Z(3)],strcat('b',linha),'Linewidth',width) %x


end

