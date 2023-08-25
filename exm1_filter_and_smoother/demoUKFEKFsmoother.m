clc
clear
close
dbstop if error
mc=1;
for i=1:mc
esterror=zeros(8,100,mc);
model=genmodel;
truth=gentruth(model);
meas=genmeas(model,truth);
est1=myUKF(model,truth,meas,1);
est2=myUKF(model,truth,meas,2);
esterror(:,:,i)=[est1.error;est1.errorsm];
end
esterror=mean(esterror,3);
x=truth.X;
sx=truth.station;
z=meas.Z;
estx=est1.Xsm;
estx2=est1.Xsm;
figure(3)
plot(x(1,:),x(3,:),'r-d','LineWidth',1)
hold on
for i=1:100
plot(z(1,i).*cos(z(2,i))+sx(1,i),z(1,i).*sin(z(2,i))+sx(3,i),'b.','LineWidth',1,'MarkerSize',3)
plot(estx(1,i),estx(3,i),'^-','Color',[0.11,0.82,0.33])
plot(estx2(1,i),estx2(3,i),'kd--')
legend('ori','meas','EKF','UKF')
grid on
pause(0.02)
end
hold off
figure(1)
plot(esterror(1,:));
hold on
plot(esterror(5,:));
legend('xerrorEKF','xerrorUKF')
hold off
figure(2)
plot(esterror(3,:));
hold on
plot(esterror(7,:));
legend('yerrorEKF','yerrorUKF')
hold on
% figure(4)
% subplot(211)
% plot(est1.error(1,:),'r--');
% hold on
% plot(est2.errorsm(1,:),'b-.');
% legend('PF','PFSM')
% subplot(212)
% plot(est1.error(3,:),'r--');
% hold on
% plot(est2.errorsm(3,:),'b-.');
% legend('PF','PFSM')