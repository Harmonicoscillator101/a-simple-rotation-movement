
% 电子在均匀 B=(0,0,Bz) 与恒定电场 E=(Ex,0,0) 下的运动（解析解）
clear; close all; clc;

% -------------------------
% 常量（SI 单位）
% -------------------------
q  = -1.602176634e-19;    % 电子电荷 (C)
m  = 9.1093837015e-31;    % 电子质量 (kg)

% -------------------------
% 磁场与电场（可修改）
% -------------------------
Bz = 1;                % 磁场强度 (T)，沿 z 方向
Ex = 5e4;                % 电场强度 (V/m)，沿 x 方向 （可为正或负）
B = [0;0;Bz];
E = [Ex;0;0];

% -------------------------
% 初始条件（位置、速度）
% -------------------------
r0 = [0; 0; 0];           % 初始位置 (m)
v0 = [1e6; 0.5e6; 0.2e6]; % 初始速度 (m/s)

% -------------------------
% 基本量
% -------------------------
if abs(Bz) < 1e-30
    error('Bz 过小或为零，无法使用本解析公式。请增大 Bz。');
end
omega = (q * Bz) / m;          % 带符号回旋角频率
omega_abs = abs(omega);

% E x B 漂移（对任意电荷均适用）
v_d = cross(E, B) / (norm(B)^2);   % 一般公式
% 对当前简单情形，v_d = [0; -Ex/Bz; 0]

% 相对速度（相对于漂移中心的初始横向速度）
v_perp0 = v0(1:2);
v_d_perp = v_d(1:2);
v_rel0 = v_perp0 - v_d_perp;   % 初始相对速度（决定回旋半径）

% Larmor 半径由相对速度决定
rL_rel = m * norm(v_rel0) / (abs(q) * abs(Bz));
Tcycl = 2*pi / omega_abs;

fprintf('Bz = %.3e T, Ex = %.3e V/m\n', Bz, Ex);
fprintf('|omega| = %.3e rad/s, 回旋周期 T = %.3e s\n', omega_abs, Tcycl);
fprintf('E×B 漂移 v_d = [%.3e, %.3e, %.3e] m/s\n', v_d);
fprintf('相对 Larmor 半径 r_L_rel ≈ %.3e m\n', rL_rel);

% -------------------------
% 时间设置
% -------------------------
nPeriods = 5;
tMax = nPeriods * Tcycl;
Nt = 2000;
t = linspace(0, tMax, Nt);

% -------------------------
% 解析解（含 E 漂移）
% 记 v_rel0 = [v_rel_x0; v_rel_y0]
% v_rel_x(t) = v_rel_x0*cos(omega t) + v_rel_y0*sin(omega t)
% v_rel_y(t) = v_rel_y0*cos(omega t) - v_rel_x0*sin(omega t)
% v_perp(t) = v_d_perp + v_rel(t)
% x(t) = x0 + v_d_x * t + (v_rel_x0/omega) * sin(omega t) - (v_rel_y0/omega)*(1 - cos(omega t))
% y(t) = y0 + v_d_y * t + (v_rel_y0/omega) * sin(omega t) + (v_rel_x0/omega)*(1 - cos(omega t))
% z(t) = z0 + v_z0 * t
% -------------------------
vx_rel0 = v_rel0(1);
vy_rel0 = v_rel0(2);
vz0 = v0(3);

% 避免 omega = 0 的奇异
if abs(omega) < 1e-30
    error('omega ≈ 0（Bz 太小），请增大 Bz。');
end

% 相对速度随时间（回旋）
vx_rel = vx_rel0 * cos(omega * t) + vy_rel0 * sin(omega * t);
vy_rel = vy_rel0 * cos(omega * t) - vx_rel0 * sin(omega * t);

% 总速度 = 漂移 + 相对速度
vx = v_d(1) + vx_rel;
vy = v_d(2) + vy_rel;
vz = vz0 * ones(size(t));

% 位置（解析积分）
x = r0(1) + v_d(1) * t + (vx_rel0/omega) .* sin(omega * t) - (vy_rel0/omega) .* (1 - cos(omega * t));
y = r0(2) + v_d(2) * t + (vy_rel0/omega) .* sin(omega * t) + (vx_rel0/omega) .* (1 - cos(omega * t));
z = r0(3) + vz0 * t;

% -------------------------
% 绘图
% -------------------------
figure('Name','电子在 E=(Ex,0,0), B=(0,0,Bz) 下的运动','Color','w','Position',[100 100 900 600]);

subplot(2,2,[1 3]);
plot3(x,y,z,'b','LineWidth',1.2); hold on;
plot3(x(1),y(1),z(1),'go','MarkerFaceColor','g');
plot3(x(end),y(end),z(end),'ro','MarkerFaceColor','r');
grid on; axis equal; view(35,25);
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('三维轨迹：回旋 + 指导心漂移');

subplot(2,2,2);
plot(x,y,'b'); hold on;
plot(x(1),y(1),'go','MarkerFaceColor','g'); plot(x(end),y(end),'ro','MarkerFaceColor','r');
axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)');
title('xy 投影（漂移使圆心沿 y 方向移动）');

subplot(2,2,4);
plot(x,z,'m'); grid on; xlabel('x (m)'); ylabel('z (m)');
title('xz 投影');


fprintf('绘图完成。注意：漂移速度大小为 |v_d| = |Ex/Bz|，方向为 -y（当 Ex>0 时）。\n');
