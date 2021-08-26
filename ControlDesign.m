clc;
%% Quasi-LPV representation
% Quasi-LPV variables
omega = sym('omega', [3 1]);
eta = sym('eta');
qq = [eta; omega]; % Time-Varying Parameters

% break the system in two subsystems
s = [];

% Rotational Subsystem
s(r).A = [(-1/2) * (skew(omega)) (1/2) * (eta * eye(3));
    zeros(3, 3) sp.Im \ skew(omega) * sp.Im];

s(r).B = [zeros(3, 3);
    inv(sp.Im)];

s(r).C = [eye(3) zeros(3, 3)];

s(r).D = zeros(3, 3);

% Translational Subsystem
s(t).A = [zeros(3, 3) eye(3);
    zeros(3, 6)];
s(t).B = [zeros(3, 3);
    sp.m \ eye(3)];

s(t).C = [eye(3) zeros(3, 3)];

s(t).D = zeros(3, 3);

%% Global System
s(g).A = blkdiag(s(r).A, s(t).A);
s(g).B = blkdiag(s(r).B, s(t).B);
s(g).C = blkdiag(s(r).C, s(t).C);
s(g).D = blkdiag(s(r).D, s(t).D);

%% Construct the internal model

% SISO synthesis
syms delta; % Time-Varying resonant frequency
Phi = []; % integrator state
Gam = []; % integrator state

% multi-resonant states
h = 3; % harmonics predicted by the control law
kr = 1; % resonant control constant

for i = 1:h
    Phi = blkdiag(Phi, [0 i * delta; -i * delta 0]);
    Gam = [Gam; kr; kr];
end

% MIMO synthesis
Phi = blkdiag(Phi, Phi, Phi); % diagonolize p times
Gam = blkdiag(Gam, Gam, Gam);

%% Polytope vertices
qq = [qq; delta]; % add IMP time varying parameter
th_max = 25 * pi / 180;
w_max = 5;
lim = [];
%              high         low
lim.eta = [1 cos(th_max / 2)];
lim.w = [w_max -w_max];
lim.delta = [0.5 0.3];

v = allcomb(lim.eta, lim.w, lim.w, lim.w, lim.delta)';

M = size(v, 2); % number of vertices

%% Augmented System
n = 6; % number of states
nu = size(Phi, 1); % number of IMP states
nn = nu + n;
sa = [];

% Decoupled models
for id = r:t
    sa(id).A = [];
    sa(id).B = [];

    A = [s(id).A zeros(n, nu);
        -Gam * s(id).C Phi];

    B = [s(id).B;
        -Gam * s(id).D];

    for i = 1:M
        sa(id).A(:, :, i) = double(subs(A, qq, v(:, i)));
        sa(id).B(:, :, i) = double(subs(B, qq, v(:, i)));
    end

end

% Unified Model
sa(g).A = [];
sa(g).B = [];

A = [blkdiag(s(1).A, s(2).A) zeros(2 * n, 2 * nu);
    -blkdiag(Gam * s(1).C, Gam * s(2).C) blkdiag(Phi, Phi)];

B = [blkdiag(s(1).B, s(2).B);
    blkdiag(-Gam * s(1).D, -Gam * s(2).D)];

for i = 1:M
    sa(g).A(:, :, i) = double(subs(A, qq, v(:, i)));
    sa(g).B(:, :, i) = double(subs(B, qq, v(:, i)));
end

%% Control Synthesis

x_bar = [sin(th_max / 2) sin(th_max / 2) sin(th_max / 2) lim.w lim.w lim.w];
p = zeros(2 * n, 2 * nn);
I = eye(2 * nn);

for k = 1:6
    p(k, :) = I(k, :) / x_bar(k);
end

K = [];
P = [];

[K.g, P.g] = lmi_with_exp_decay(sa(g).A, sa(g).B, p, 0.05, 15);

