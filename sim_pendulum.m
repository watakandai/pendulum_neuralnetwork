function f = sim_pendulum(pop)

global M m L J Dc Dth g dt
global hid_num out_num in_num IH_num
global pops

%% GAの遺伝子情報をNNの重みとして利用
X = [-1 1;-3.14 3.14];
net = newff(minmax(X),[hid_num,out_num],{'tansig','hardlim'});
W_IH = zeros(hid_num,in_num);
W_HO = zeros(out_num,hid_num);
for ii = 1:pops
    for n1=1:hid_num
        n2 = (n1-1)*in_num + 1;
        for n3 = 1:in_num
            n4 = n2 + n3 - 1;
            W_IH(n1,n3) = pop(ii,n4);
        end
    end
    
    for n1 = 1:out_num
        n2 = (n1-1)*hid_num + 1;
        for n3 = 1:hid_num
            n4 = n2 + n3 - 1;
            W_HO(n1,n3) = pop(ii,n4+IH_num);
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
    
    states = {'z' 'th' 'dz' 'dth'};
    inputs = {'u'};
    outputs = {'z'; 'th'};
    sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);
    
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
        u = -K*X;
        
        % Data Strage
        X_data(:,i) = X;
%         U_data(1,i) = U;
        u_data(1,i) = u;
        time(i,1) = dt * i;
        
        % Runge-Kutta
        % Linear
%         X1 = X;        k1 = dt*(A*X1+B*u);
%         X2 = X+k1/2;   k2 = dt*(A*X2+B*u);
%         X3 = X+k2/2;   k3 = dt*(A*X3+B*u);
%         X4 = X+k3;     k4 = dt*(A*X4+B*u);
%         % ----- (i+1) step ----- %
%         X = X+(k1+2*k2+2*k3+k4)/6;

        % Nonlinear
        dX1 = getdX(X, u)*dt;
        dX2 = getdX(X+dX1/2, u)*dt;
        dX3 = getdX(X+dX2/2, u)*dt;
        dX4 = getdX(X+dX3, u)*dt;
        X = X+(dX1+2*dX2+2*dX3+dX4)/6;

        if abs(X(1,1)) > 1.0
            break;
        end
        %% animation
%         cla;
%         figure(2)
%         line([-1 1], [0 0]);
%         line([X(1,1), X(1,1)+2*L*sin(X(3,1))], [0 2*L*cos(X(3,1))]);
%         axis([-1 1 -0.5 0.5]); xlabel('x [m]'); ylabel('z [m]'); hold on;
%         rectangle('Position',[X(1,1)-0.05 -0.025 0.1 0.05],'Curvature',0.2); grid on;
%         drawnow;
    end
    f(ii,1) = sum(cos(X(2,:)));
end