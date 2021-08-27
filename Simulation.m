%% Creates Object
Ts = 1e-3;

sp = StewartPlatform;
InitialCondtions = [];
InitialCondtions.q = Quaternion.ea2q([0 0 0]);
InitialCondtions.w = [0; 0; 0];
InitialCondtions.p = sp.tT + 0;
InitialCondtions.v = [0; 0; 0];

for freq = [0.3  0.5]

    sp = StewartPlatform(Ts, InitialCondtions);
    %% Tests
    x = [sp.x(1).q.e; sp.x(1).w; sp.x(1).p; sp.x(1).v];

    InitialState = [sp.x(1).q.n;
                sp.x(1).w;
                freq];

    Sw = double(subs(blkdiag(Phi, Phi), qq, InitialState)); % B \Omega

    M = double(subs([s(g).A + s(g).B * K.g(:, 1:12) s(g).B * K.g(:, 13:end);
                blkdiag(Gam, Gam) * s(g).C + blkdiag(Gam, Gam) * s(g).D * K.g(:, 1:12) blkdiag(Phi, Phi) + blkdiag(Gam, Gam) * s(g).D * K.g(:, 13:end)], qq, InitialState));
    U = [];

    for i = 1:h
        U = [U 1 0];
    end

    U = blkdiag(U, U, U, U, U, U);

    Bw = s(g).B * U;

    Dw = zeros(6, 12 * h);

    Omega = [Bw; blkdiag(Gam, Gam) * Dw];

    Pi_Sigma = sylvester(M, -Sw, Omega);

    Sigma = Pi_Sigma(13:end, :);

    w = [5.6e-4 * ones(6 * h, 1); ones(6 * h, 1)];
    z = [x; -Sigma * w];

    [z' * P.g * z]

    if [z' * P.g * z] > 1
        disp 'Invalid intial values';
        return
    end

    disp 'Valid initial values';
    disp 'Starting simulation...';

    %% Simulation
    Ttotal = 60;

    Ref = []; Ref.r = [0; 0; 0]; Ref.t = sp.tT;

    eg = []; eg.r = LinearSystem(Phi, [], U(1:3, 1:6 * h), [], sp.Ts, w(1:6 * h)); eg.t = LinearSystem(Phi, [], U(4:6, 6 * h + 1:end), [], sp.Ts, w(6 * h + 1:end));

    xi = []; xi.r = LinearSystem (Phi, Gam, eye(size(Phi, 1)), [], sp.Ts); xi.t = LinearSystem (Phi, Gam, eye(size(Phi, 1)), [], sp.Ts);

    x = [];

    x.r = [sp.x(1).q.e * sign(sp.x(1).q.n); sp.x(1).w];

    x.t = [sp.x(1).p; sp.x(1).v];

    sp = StewartPlatform(Ts, InitialCondtions);

    for i = 1:round(Ttotal / sp.Ts)

        % rotational subsystem - s0
        x.r = [sp.x(i).q.e; sp.x(i).w];
        e.r = -s(r).C * x.r;

        % translational subsystem - s1
        x.t = [sp.x(i).p - Ref.t; sp.x(i).v]; % Ref.t is a offset, we changed the coordinates to the origin be at home position
        e.t = -s(t).C * x.t;

        % control law
        u = K.g(:, 1:12) * [x.r; x.t] + K.g(:, 13:end) * [xi.r.Sys_Dynamics(e.r, delta, freq); xi.t.Sys_Dynamics(e.t, delta, freq)];

        D = [eg.r.Sys_Dynamics([], delta, freq);
                            eg.t.Sys_Dynamics([], delta, freq)];

        % Simulation
        sp.Sys_Dynamics(sp.J' \ u, sp.J' \ D);

    end

    %%
    figure
    subplot(3, 2, 2)
    p = plot(sp.t, [sp.x.p] - Ref.t, 'linewidth', 1.5);
    title('position error')
    grid on;
    xlabel('t(s)');
    ylabel('p(m)')
    xy_axis = axis;
    axis([0 Ttotal xy_axis(3) xy_axis(4)]);
    set(p(2), 'LineStyle', '-.');
    set(p(3), 'LineStyle', '--');

    q = [sp.x.q];
    subplot(3, 2, 1)
    p = plot(sp.t, [q.e], 'linewidth', 1.5);
    title('orientation error')
    grid on;
    xlabel('t(s)');
    ylabel('\epsilon')
    xy_axis = axis;
    axis([0 Ttotal xy_axis(3) xy_axis(4)]);
    set(p(2), 'LineStyle', '-.');
    set(p(3), 'LineStyle', '--');

    subplot(3, 2, 3)
    p = plot(sp.t, sp.u(1:3, :), 'linewidth', 1.5);
    hold on
    grid on;
    xlabel('t(s)');
    ylabel('u_\tau (Nm)')
    title('torque input')
    xy_axis = axis;
    axis([0 Ttotal xy_axis(3) xy_axis(4)]);
    set(p(2), 'LineStyle', '-.');
    set(p(3), 'LineStyle', '--');

    subplot(3, 2, 4)
    p = plot(sp.t, sp.u(4:6, :), 'linewidth', 1.5);
    hold on
    grid on;
    xlabel('t(s)');
    ylabel('u_F (N)')
    title('force input')
    xy_axis = axis;
    axis([0 Ttotal xy_axis(3) xy_axis(4)]);
    set(p(2), 'LineStyle', '-.');
    set(p(3), 'LineStyle', '--');

    subplot(3, 2, 5)
    p = plot(sp.t, sp.d(1:3, :), 'linewidth', 1.5);
    hold on;
    grid on;
    xlabel('t(s)');
    ylabel('d_\tau (Nm)')
    title('pertubation torque')
    xy_axis = axis;
    axis([0 Ttotal xy_axis(3) xy_axis(4)]);
    set(p(2), 'LineStyle', '-.');
    set(p(3), 'LineStyle', '--');

    subplot(3, 2, 6);
    p = plot(sp.t, sp.d(4:6, :), 'linewidth', 1.5);
    hold on;
    grid on;
    xlabel('t(s)');
    ylabel('d_F (N)');
    title('perturbation force');
    xy_axis = axis;
    axis([0 Ttotal xy_axis(3) xy_axis(4)]);
    set(p(2), 'LineStyle', '-.');
    set(p(3), 'LineStyle', '--');

    l = legend('x-axis', 'y-axis', 'z-axis');
    set(l, ...
        'Position', [0.48267447598653 0.607536196568469 0.0585651531008292 0.0721966186609869]);

end
