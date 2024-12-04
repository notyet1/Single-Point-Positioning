data = readtable('XYZerror.csv');  

x_error = data.X_error;  % X方向误差
y_error = data.Y_error;  % Y方向误差
z_error = data.Z_error;  % Z方向误差

n = length(x_error);  

time = 1:n;  

figure;  

% 绘制 X 方向误差
subplot(3, 1, 1);  
plot(time, x_error, 'r');  
title('X Direction Error');
xlabel('Iterations');
ylabel('X Error (m)');
grid on;

% 绘制 Y 方向误差
subplot(3, 1, 2);  
plot(time, y_error, 'g');  
title('Y Direction Error');
xlabel('Iterations');
ylabel('Y Error (m)');
grid on;

% 绘制 Z 方向误差
subplot(3, 1, 3);  
plot(time, z_error, 'b');  
title('Z Direction Error');
xlabel('Iterations');
ylabel('Z Error (m)');


% 测站的地理坐标（单位为度）
%纬度: 30.52781930°
%经度: 114.35656423°
%高度: 80.804 m
latitude = 30.5;    
longitude = 114.3; 

% 将经纬度转换为弧度
phi = deg2rad(latitude);
lambda = deg2rad(longitude);

% 提取 XYZ 误差数据
x_error = data.X_error;
y_error = data.Y_error;
z_error = data.Z_error;

% 初始化 NEU 误差向量
n_error = zeros(length(x_error), 1);
e_error = zeros(length(x_error), 1);
u_error = zeros(length(x_error), 1);

% 计算 NEU 误差
for i = 1:length(x_error)
    % ECEF 误差向量
    dXYZ = [x_error(i); y_error(i); z_error(i)];
    % 转换矩阵
    T = [-sin(phi) * cos(lambda), -sin(phi) * sin(lambda), cos(phi);
         -sin(lambda),             cos(lambda),            0;
          cos(phi) * cos(lambda),  cos(phi) * sin(lambda), sin(phi)];
    % 计算 NEU 误差
    dNEU = T * dXYZ;
    % 存储 NEU 误差
    n_error(i) = dNEU(1);  
    e_error(i) = dNEU(2);  
    u_error(i) = dNEU(3);  
end


time = 1:length(n_error);

figure;

% 绘制 N 方向误差
subplot(3, 1, 1);
plot(time, n_error, 'r');
title('North Direction Error');
xlabel('Iterations');
ylabel('N Error ');
grid on;

% 绘制 E 方向误差
subplot(3, 1, 2);
plot(time, e_error, 'g');
title('East Direction Error');
xlabel('Iterations');
ylabel('E Error ');
grid on;

% 绘制 U 方向误差
subplot(3, 1, 3);
plot(time, u_error, 'b');
title('Up Direction Error');
xlabel('Iterations');
ylabel('U Error ');
grid on;

% 调整子图布局
tight_layout();
