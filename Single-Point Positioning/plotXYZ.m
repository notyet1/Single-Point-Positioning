data = readtable('XYZerror.csv');  

x_error = data.X_error;  % X�������
y_error = data.Y_error;  % Y�������
z_error = data.Z_error;  % Z�������

n = length(x_error);  

time = 1:n;  

figure;  

% ���� X �������
subplot(3, 1, 1);  
plot(time, x_error, 'r');  
title('X Direction Error');
xlabel('Iterations');
ylabel('X Error (m)');
grid on;

% ���� Y �������
subplot(3, 1, 2);  
plot(time, y_error, 'g');  
title('Y Direction Error');
xlabel('Iterations');
ylabel('Y Error (m)');
grid on;

% ���� Z �������
subplot(3, 1, 3);  
plot(time, z_error, 'b');  
title('Z Direction Error');
xlabel('Iterations');
ylabel('Z Error (m)');


% ��վ�ĵ������꣨��λΪ�ȣ�
%γ��: 30.52781930��
%����: 114.35656423��
%�߶�: 80.804 m
latitude = 30.5;    
longitude = 114.3; 

% ����γ��ת��Ϊ����
phi = deg2rad(latitude);
lambda = deg2rad(longitude);

% ��ȡ XYZ �������
x_error = data.X_error;
y_error = data.Y_error;
z_error = data.Z_error;

% ��ʼ�� NEU �������
n_error = zeros(length(x_error), 1);
e_error = zeros(length(x_error), 1);
u_error = zeros(length(x_error), 1);

% ���� NEU ���
for i = 1:length(x_error)
    % ECEF �������
    dXYZ = [x_error(i); y_error(i); z_error(i)];
    % ת������
    T = [-sin(phi) * cos(lambda), -sin(phi) * sin(lambda), cos(phi);
         -sin(lambda),             cos(lambda),            0;
          cos(phi) * cos(lambda),  cos(phi) * sin(lambda), sin(phi)];
    % ���� NEU ���
    dNEU = T * dXYZ;
    % �洢 NEU ���
    n_error(i) = dNEU(1);  
    e_error(i) = dNEU(2);  
    u_error(i) = dNEU(3);  
end


time = 1:length(n_error);

figure;

% ���� N �������
subplot(3, 1, 1);
plot(time, n_error, 'r');
title('North Direction Error');
xlabel('Iterations');
ylabel('N Error ');
grid on;

% ���� E �������
subplot(3, 1, 2);
plot(time, e_error, 'g');
title('East Direction Error');
xlabel('Iterations');
ylabel('E Error ');
grid on;

% ���� U �������
subplot(3, 1, 3);
plot(time, u_error, 'b');
title('Up Direction Error');
xlabel('Iterations');
ylabel('U Error ');
grid on;

% ������ͼ����
tight_layout();
