% demo for smooth
clc
clear
close
%% 1.initiation
t=100;
T=3 ; 
Num=round(t/T);
X=[100;1;100;1];
F=[1 T 0 0;
   0 1 0 0;
   0 0 1 T;
   0 0 0 1];
turn=pi/60;
F2=[1,0,(sin(-turn*T)/-turn),(-(1-cos(-turn*T))/-turn);
      0,1,(1-cos(-turn*T))/-turn,(sin(-turn*T)/-turn);
      0,0,cos(-turn*T),(-sin(-turn*T));
      0,0,sin(-turn*T),cos(-turn*T)];
H=[1 0 0 0;
   0 0 1 0];
% Hn=@(x)[norm(x([1,3]));atan2(X(3),X(1))];
P=eye(4)*5;
Q=diag([0.5;1;0.5;1]).^2;
R=diag([20;20]).^2;
Xe=[50;2;50;1];Xecorr=Xe;
Xestore=repmat(zeros(size(Xe)),[1,Num]);
GT=Xestore;
Pstore=repmat(zeros(size(P)),[1,1,Num]);
PPstrore=repmat(zeros(size(P)),[1,1,Num]);
xsmooth=Xestore;Psmooth=Pstore;PPsmooth=Psmooth;
%% 2.Forward filtering
for i=1:Num
    X=F2*X+sqrt(Q)*randn(4,1)+[1 0;0.1,0.1;0 1;0.1,0.1]*[0.2;0.3];
    Z=H*X+sqrt(R)*randn(2,1);
 [Xe,~,P,PP]=myKF('X',Xe,'P',P,'Z',Z,'Q',Q,'R',R,'F',F);
 Pstore(:,:,i)=P;
 Xestore(:,i)=Xe;
 PPstrore(:,:,i)=PP;
 GT(:,i)=X;
end
%% Backward smoothing
xsmooth(:,Num) = Xestore(:,Num);
Psmooth(:,:,Num) = Pstore(:,:,Num);
for i=Num-1:-1:1
   [xsmooth(:,i), Psmooth(:,:,i), PPsmooth(:,:,i+1)]=myKFsmoother(xsmooth(:,i+1), Psmooth(:,:,i+1), ...
    Xestore(:,i), Pstore(:,:,i),  Pstore(:,:,i+1), PPstrore(:,:,i+1), F, Q);
end
errorKF=sqrt((GT-Xestore).^2);
errorKFsmooth=sqrt((GT-xsmooth).^2);
figure(1)
plot(GT(1,:),GT(3,:),'k-','LineWidth',1.5);
hold on
scatter(GT(1,1),GT(3,1),'k');scatter(GT(1,end),GT(3,end),'cyan');
plot(Xestore(1,:),Xestore(3,:),'-','LineWidth',1.5);
plot(xsmooth(1,:),xsmooth(3,:),'-','LineWidth',1.5);
legend('GT','start','end','KF','KFsmoother','Location','best')
title('跟踪效果')
axis 'auto xy'
figure(2)
subplot(211)
plot(1:Num,errorKFsmooth(1,:),'r','LineWidth',1.5)
hold on
plot(1:Num,errorKF(1,:),'b','LineWidth',1.5)
title('x坐标误差')
xlabel('m','Interpreter','latex')
ylabel('m','Interpreter','latex')
legend('kfsmooth','kf')
subplot(212)
plot(1:Num,errorKFsmooth(3,:),'r','LineWidth',1.5)
hold on
plot(1:Num,errorKF(3,:),'b','LineWidth',1.5)
title('y坐标误差')
xlabel('m','Interpreter','latex')
ylabel('m','Interpreter','latex')
legend('kfsmooth','kf')