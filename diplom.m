%% 
clc; close all; clear;
format compact
format shortG


% Глобальные переменные для ПИД-регуляторов
global error_u error_r u_dPrev r_dPrev time_prev;
error_u = 0;
error_r = 0;
u_dPrev = 0;
r_dPrev = 0;
time_prev = 0;

% Координаты точки старта (м)
x0 = 0;
y0 = 0;
% Координаты цели (м)
xf = 20;
yf = 20;

% Начальный угол к цели
phi0 = atan2((yf - y0), (xf - x0)) + 0.08;
rad2deg(phi0)

% Начальная скорость (м/с)
V0 = 1.0; 
% Начальные условия: [x, y, phi, u_r, v_r, r]
F0 = [x0 y0 phi0 V0 0 0];
% Конечная точка
Pf = [xf yf 0 0 0 0];

% Настройки ПИД-регулятора для tau_u_r
Kp = 200;
Ki = 30;
Kd = 0; 
Kpid_u = [Kp Ki Kd];
% Настройки ПИД-регулятора для tau_phi_n
Kp = 100;
Ki = 0;
Kd = 0; 
Kpid_r = [Kp Ki Kd];

% Скорость течения
VC = 0.05;
beta_c = 0.1;
rad2deg(beta_c)
VC/V0

% Расстояние до цели - начальное (м)
Do = sqrt((xf-x0)^2 + (yf-y0)^2);
% Расстояние для снижения скорости (м)
Dk = 0.5; 
% Коэффициент снижения скорости
a = V0/Dk;
% Степень в уравнении снижения скорости
s = 1;
Vs = VC;

% ======================
% Метод Эйлера
% ======================
dt = 0.01;
tspan = 0:dt:30;
n = length(tspan);
F = zeros(n, 6);
F(1, :) = F0;
t = tspan;

% Коэффициенты динамики
m11 = 47.5;
m22 = 94.1;
m33 = 13.6;
m23 = 5.2;
m32 = 5.2;
d11 = 13.5;
d22 = 50.2;
d33 = 27.2;
d23 = 41.4;
d32 = 17.3;

% Матрицы
M = [m11, 0, 0; 0, m22, m23; 0, m32, m33];
D = [d11, 0, 0; 0, d22, d23; 0, d32, d33];
invM = M^-1;

% Цикл метода Эйлера
for i = 1:n-1
    % Текущее состояние
    x_n = F(i, 1);
    y_n = F(i, 2);
    phi_n = F(i, 3);
    u_r = F(i, 4);
    v_r = F(i, 5);
    r = F(i, 6);

    % Расстояние до цели
    D = sqrt((xf-x_n)^2 + (yf-y_n)^2);
    V_aim = V0;
    if D < Dk
        V_aim = (Dk/Dk^s)*a*D^s;
    end

    % Скорость течения
    u_c = VC * cos(beta_c - phi_n);
    v_c = VC * sin(beta_c - phi_n);
    Vc = [u_c; v_c; 0];
    Vr = [u_r; v_r; r];
    V = Vr + Vc;

    % Матрицы уравнений
    C = [0, 0, -m22*v_r-m23*r;
         0, 0, m11*u_r;
         m22*v_r+m23*r, -m11*u_r, 0];
    R = [cos(phi_n), -sin(phi_n), 0;
         sin(phi_n), cos(phi_n), 0;
         0, 0, 1];

    % Управление скоростью u_r
    dVrdt = -invM*(C*Vr + D*Vr);
    e_u = V_aim - u_r;
    if t(i) == time_prev
        u_dDot = 0;
    else
        u_dDot = (u_r - u_dPrev) / (t(i) - time_prev);
    end
    u_dPrev = u_r;
    tau_u = Kpid_u(1)*e_u + Kpid_u(3)*(u_dDot - dVrdt(1)) + Kpid_u(2)*error_u;
    error_u = error_u + e_u*(t(i) - time_prev);
    error_u = min(max(error_u, -100), 100);

    % Управление углом phi_n
    phi_aim = atan2((yf - y_n), (xf - x_n));
    e_r = phi_aim - phi_n;
    if t(i) == time_prev
        r_dDot = 0;
    else
        r_dDot = (phi_n - r_dPrev) / (t(i) - time_prev);
    end
    r_dPrev = phi_n;
    tau_r = Kpid_r(1)*e_r + Kpid_r(3)*(r_dDot - r) + Kpid_r(2)*error_r;
    error_r = error_r + e_r*(t(i) - time_prev);
    error_r = min(max(error_r, -100), 100);
    % ==========
    % Ограничение всплексков управляющего воздейсвия
    % до реальных значений
    % ==========
    tau_u_max = 200;
    if abs (tau_u) > tau_u_max
        tau_u = sign (tau_u) * tau_u_max;
    end
 

    tau_r_max = 50;
    if abs (tau_r) > tau_r_max
        tau_r = sign (tau_r) * tau_r_max;
    end


    % Управляющие воздействия
    TAU = [tau_u; 0; tau_r];

    % Производные
    dETA_n = R*V;
    dVr = invM*(TAU - C*Vr - D*Vr);
    dFdt = [dETA_n(1); dETA_n(2); dETA_n(3); dVr(1); dVr(2); dVr(3)];

    % Проверка на корректность
    if any(isnan(dFdt)) || any(isinf(dFdt))
        fprintf('Ошибка: NaN или Inf в dFdt на шаге %d, t=%.2f\n', i, t(i));
        break;
    end

    % Обновление состояния
    F(i+1, :) = F(i, :) + dt * dFdt';
    if any(isnan(F(i+1, :))) || any(isinf(F(i+1, :)))
        fprintf('Ошибка: NaN или Inf в F на шаге %d, t=%.2f\n', i, t(i));
        break;
    end

    time_prev = t(i);
end

% Проверка данных
if any(isnan(F(:))) || any(isinf(F(:)))
    error('Ошибка: Массив F содержит NaN или Inf. Графики не могут быть построены.');
else
    disp('Данные в F корректны, строим графики.');
end

% Извлечение результатов
x = F(:,1);
y = F(:,2);
phi_n = F(:,3);
u_r = F(:,4);
v_r = F(:,5);
r = F(:,6);

% Скорость течения
u_c = VC * cos(beta_c - phi_n);
v_c = VC * sin(beta_c - phi_n);
Vc = [u_c, v_c, zeros(length(u_r),1)];
Vr = [u_r, v_r, r];
V = Vr + Vc;
u = V(:,1);
v = V(:,2);

% Угол направления к цели
phi_aim = atan2((yf - y), (xf - x));

% Ошибки
Ex = xf - x;
Ey = yf - y;
Ephi = atan2((yf - y), (xf - x)) - phi_n;
ED = sqrt(Ex.^2 + Ey.^2);

% Графики скорости ПА во времени
figure('Name', 'Скорости AUV');
m = 2; n = 2;
subplot(m,n,1); plot(t, u_r, 'k', t, u, 'r'); grid on;
xlabel('t'); ylabel('u_r, u'); title('Скорость AUV: u_r, u');
lgd = legend('u_r','u');
lgd.Location = 'best';

subplot(m,n,2); plot(t, v_r, 'k', t, v, 'r'); grid on;
xlabel('t'); ylabel('v_r, v'); title('Скорость AUV: v_r, v');
lgd = legend('v_r','v');
lgd.Location = 'best';

subplot(m,n,3); plot(t, r, 'k'); grid on;
xlabel('t'); ylabel('r'); title('Угловая скорость r');

% Графики параметров ПА во времени
figure('Name', 'Параметры AUV');
m = 2; n = 2;
subplot(m,n,1); plot(t, x, 'k'); grid on;
xlabel('t'); ylabel('x'); title('Координата x');

subplot(m,n,2); plot(t, y, 'k'); grid on;
xlabel('t'); ylabel('y'); title('Координата y');

subplot(m,n,3); plot(t, phi_n, 'k', t, phi_aim, 'r'); grid on;
xlabel('t'); ylabel('\phi'); title({'Курсовой угол и угол до цели'});
lgd = legend('курс','угол от СТЗ');
lgd.Location = 'best';

subplot(m,n,4); plot(t, r, 'k-'); grid on;
xlabel('t'); ylabel('r'); title('Скорость r');

% Графики ошибок во времени
figure('Name', 'Ошибки AUV');
m = 2; n = 2;
subplot(m,n,1); plot(t, Ex, 'k'); grid on;
xlabel('t'); ylabel('E_x'); title('Ошибка по x');

subplot(m,n,2); plot(t, Ey, 'k'); grid on;
xlabel('t'); ylabel('E_y'); title('Ошибка по y');

subplot(m,n,3); plot(t, Ephi, 'k'); grid on;
xlabel('t'); ylabel('E_\phi'); title('Ошибка по углу');

subplot(m,n,4); plot(t, ED, 'k.'); grid on;
xlabel('t'); ylabel('D'); title('Расстояние до цели');
axis([0 inf -1 15]);

% График движения ПА
figure('Name', 'Траектория AUV');
l = length(tspan)/Do;
plot(x, y, 'k-'); hold on;
plot(F(1:l:end,1), F(1:l:end,2), 'k.');
plot(F0(:,1), F0(:,2), 'k*', 'LineWidth', 2, 'MarkerSize', 10); grid on;
plot(Pf(:,1), Pf(:,2), 'rs', 'LineWidth', 2, 'MarkerSize', 15); grid on;
lgd = legend('','старт','цель');
title(lgd, 'Обозначения:');
lgd.Location = 'best';
xlabel('x'); ylabel('y');
title('График движения ПА в координатных осях');
axis equal;