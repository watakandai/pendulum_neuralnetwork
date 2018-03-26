% genetic algorithm
% 81722092
% 渡邊幹大

clear
close all

%% グローバル変数の定義
global bound rng
global hid_num out_num in_num IH_num
global pops

%% NN パラメータ
% 
%  PREFORMATTED
%  TEXT
% 
in_num  = 2;                        % 入力層のニューロン数
hid_num = 5;                        % 隠れ層のニューロン数
out_num = 1;                        % 出力層のニューロン数
IH_num  = in_num  * hid_num;        % 入力層から隠れ層への重み
HO_num  = hid_num * out_num;        % 隠れ層から出力層への重み
IHO_num = IH_num + HO_num;          % 全重み

I_in  = zeros(in_num,1);            % 入力層への入力
I_hid = zeros(hid_num,1);           % 入力層からの出力、隠れ層への入力
I_out = zeros(out_num,1);           % 隠れ層からの出力、出力層への入力
%% パラメータの初期化
pops=30;                      % 個体数
maxgen=150;                    % 世代数
crossp=0.8;                   % 交叉確率
mutatep=0.35;                 % 突然変異確率
bound=10*ones(IHO_num, 2); bound(:,1)=-10;
% bound = [-5 -4; -5 -4; 0 1; -1 0;
%          -4 -3; -1 0; 1 2; -1 0;
%          -6 -5; -7 -6; -1 0; 0 1;
%          -3 -2; -5 -4; 1 2; 0 1;
%          0 1; -7 -6; -2 -1; 0 1;
%          0 1; -1 1; 0 1; -1 0; 0 1];
numvar=size(bound,1);         % 染色体の長さ
rng=(bound(:,2)-bound(:,1))'; % 変数の範囲

pop=zeros(pops,numvar);       % 個体の初期化
pop(:,1:numvar)=(ones(pops,1)*rng).*(rand(pops,numvar))+(ones(pops,1)*bound(:,1)'); % 個体の生成
   
%% 世代の開始
for it=1:maxgen
    it
    fpop=sim_pendulum(pop);    % 適応度の計算
%     fpop=multipeak(pop);
    [cs,inds]=max(fpop);    % エリート　cs:最大値  inds:順番
    bchrom=pop(inds,:);     % エリートの値の格納
    
    % 選択
    toursize=5;
    players=ceil(pops*rand(pops,toursize)); % 適応度の組み合わせ
    scores=fpop(players);
    [a,m]=max(scores');
    pind=zeros(1,pops);
    for i=1:pops
        pind(i)=players(i,m(i));	
        parent(i,:)=pop(pind(i),:);	
    end
    
    % 交叉	
    child=cross(parent,crossp);
    
    % 突然変異
    pop=mutate(child,mutatep);
     
%     mm=sim_pendlum(pop);	
%     maxf(it)=max(mm);	
%     meanf(it)=mean(mm);
%     
%     [bfit,bind]=max(mm);
%     bsol=pop(bind,:);
%     
    % 図の作成
%     [x,y]=meshgrid([-1:0.05:1]);
%     r=sqrt(x.^2+y.^2);
%     s=sqrt((x-0.5).^2+y.^2);
%     ss=sqrt((x-0.8).^2+y.^2);
%     fff=exp(-2*r.^2)+2*exp(-1000*s.^2)+3*exp(-1000*ss.^2);
%     cla
%     mesh(x,y,fff),hold on

    % 各個体の値をプロット
%     plot3( pop(:,1),pop(:,2),mm,'r+');
    
%     % 適応度の高い個体の値をプロット
%     plot3( bsol(1),bsol(2),bfit,'md');
%     axis([-1.5 ,1.5,-1.5 ,1.5])
%     xlabel(bsol(1))
%     ylabel(bsol(2))
%     zlabel(bfit)
%     title(['Generation=',num2str(it)])
%     pause(0)

%     pop(inds,:)=bchrom;
    pop(1,:)=bchrom;
    save pop.mat pop
end
    
disp(['x=',num2str(bsol(1))])       %  NUM2STR   数値を文字列に変換
disp(['y=',num2str(bsol(2))])
disp(['z=',num2str(bfit)])