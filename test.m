clear
close all
global M m L J Dc Dth g

load pop.mat

in_num  = 2;                        % 入力層のニューロン数
hid_num = 5;                        % 隠れ層のニューロン数
out_num = 1;                        % 出力層のニューロン数
IH_num  = in_num  * hid_num;        % 入力層から隠れ層への重み
HO_num  = hid_num * out_num;        % 隠れ層から出力層への重み
IHO_num = IH_num + HO_num;          % 全重み


X = [-1 1;-3.14 3.14];
net = newff(minmax(X),[hid_num,out_num],{'tansig','hardlim'});
W_IH = zeros(hid_num,in_num);
W_HO = zeros(out_num,hid_num);
for n1=1:hid_num
    n2 = (n1-1)*in_num + 1;
    for n3 = 1:in_num
        n4 = n2 + n3 - 1;
        W_IH(n1,n3) = pop(1,n4);
    end
end

for n1 = 1:out_num
    n2 = (n1-1)*hid_num + 1;
    for n3 = 1:hid_num
        n4 = n2 + n3 - 1;
        W_HO(n1,n3) = pop(1,n4+IH_num);
    end
end
net.IW{1,1} = W_IH; % popの値を入れる
net.LW{2,1} = W_HO; % popの値を入れる
net.b{1} = zeros(hid_num,1); % popの値を入れる
net.b{2} = zeros(out_num,1); % popの値を入れる

%% Parameter (Model)
M   = 0.50;     % kg
m   = 0.10;     % kg
L   = 0.20;     % m
J   = 1.3e-3;   % kgm^2
Dc  = 1.0e-6;   % kg/s
Dth = 1.0e-6;   % kgm^2/s
g   = 9.8;      % m/s^2

%% Parameter (Simulation)
dt = 0.005;
sim_start = 0;
sim_stop = 10;
sim_num = (sim_stop-sim_start)/dt;
t = 0:dt:sim_stop-dt;     % 0 ~ 10 sec;
%% Linearlized model
Mmat = [ M+m  m*L; m*L  J+m*L^2];
Cmat = [ Dc 0; 0 Dth];
Gmat = [ 0 0; 0 -m*L*g];
Umat = [ 1; 0];

% dX = A*X + B*u, X = [z th dz dth]
A = [ zeros(2)    eye(2)
    -Mmat\Gmat  -Mmat\Cmat];
B = [ zeros(2,1) ; Mmat\Umat  ];

C = [1 0 0 0;
    0 1 0 0];
D = [0;
    0];
%% Initialization
X = [0 pi 0 0]';      % X = [z th dz dth]
u = 0;
X_data = zeros(4,sim_num);
u_data = zeros(1,sim_num);

%% 最適レギュレータ法による状態フィードバックゲインの設計
Q=diag([10,1,1,1]);
r=10;
K=lqr(A,B,Q,r);
eig(A-B*K);

%% 数値シミュレーション
for i=1:sim_num
    n = sim(net, X(1:2,1));
    % NN Controller
    switch n
        case 1
            Q =diag([10 10000 1 30]);
            r =1;
            K =lqr(A,B,Q,r);
        case 0
            if X(2,1)-pi >= 0
                K = [1 1 1 1];
            else
                K = -[1 1 1 1];
            end
    end
    u=-K*X;
    
    % Data Strage
    X_data(:,i) = X;
    u_data(1,i) = u;
    time(i,1) = dt * i;
    
    % Nonlinear
    dX1 = getdX(X, u)*dt;
    dX2 = getdX(X+dX1/2, u)*dt;
    dX3 = getdX(X+dX2/2, u)*dt;
    dX4 = getdX(X+dX3, u)*dt;
    X = X+(dX1+2*dX2+2*dX3+dX4)/6;
    
    if abs(X(1,1)) > 1.0
%         break;
    end
end

%%
figure(2)
subplot(511)
plot(t,X_data(1,:)); grid on;
xlabel('time[s]'); ylabel('z [m]');

subplot(512)
plot(t,X_data(2,:)); grid on;
xlabel('time[s]'); ylabel('\theta [rad]');

subplot(513)
plot(t,X_data(3,:)); grid on;
xlabel('time[s]'); ylabel('dz [m/s]');

subplot(514)
plot(t,X_data(4,:)); grid on;
xlabel('time[s]'); ylabel('d\theta [rad/s]');

subplot(515)
plot(t,u_data(1,:)); grid on;
xlabel('time[s]'); ylabel('u [N]');

legend('{Linear}');
%% animation
currFrame(ceil(length(t))) = struct('cdata',[],'colormap',[]);

figure(1)
for i=1:length(t)
    cla
    X = X_data(:,i);
    line([-1 1], [0 0]); grid on; hold on;
    line([X(1,1) X(1,1)+2*L*sin(X(2,1))], [0 2*L*cos(X(2,1))]); hold on;
    rectangle('Position', [X(1,1)-0.05 -0.025 0.1 0.05])
    axis([-1 1 -0.5 0.5]); xlabel('z[m]'); ylabel('y[m]');
    
    currFrame(i) = getframe(gcf);
    drawnow;
end

vidObj = VideoWriter('sim_pendulum');
open(vidObj)
for i=1:length(t)
    writeVideo( vidObj, currFrame(i) );
end
close( vidObj );

