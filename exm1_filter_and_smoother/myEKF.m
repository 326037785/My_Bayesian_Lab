function est=myEKF(model,truth,meas)
k=model.K;F=model.F1;
x=truth.X;R=model.R;
a=max(x(1,1:4))+10;b=max(x(3,1:4))+10;
est.X=zeros(model.xdim,k);
est.error=zeros(model.xdim,k);
Q=model.Q;
xstart=[randi([round(a-a/16),round(a+a/16)],1),20,randi([round(b-b/16),round(b+b/16)],1),20]';
P=diag([100,20,100,20])^2;
I=eye(model.xdim);
for i=1:k
    xstart=F*xstart;
    Pre=F*P*F'+Q;
    H=ekfhfun(xstart([1,3]),truth.station([1,3],i));
    K=Pre*H'*inv(H*Pre*H'+R);
    xstart=xstart+K*(meas.Z(:,i)-hfun(truth.station([1,3],i),xstart([1,3])));
    PP=(I-K*H)*F*P;
    P=(I-K*H)*Pre;
    est.P(:,:,i)=P;
    est.PP(:,:,i)=PP;
    est.X(:,i)=xstart;
    est.error(:,i)=sqrt((xstart-x(:,i)).^2);
end
est.xsmooth=zeros(size(est.X));
est.Psmooth=zeros(size(est.P));
est.PPsmooth=est.Psmooth;
est.xsmooth(:,k) = est.X(:,k);
est.Psmooth(:,:,k) = est.P(:,:,k);
est.errorsm(:,k)=est.error(:,k);
for i=k-1:-1:1
   [est.xsmooth(:,i), est.Psmooth(:,:,i), est.PPsmooth(:,:,i+1)]=myKFsmoother(est.xsmooth(:,i+1), est.Psmooth(:,:,i+1), ...
    est.X(:,i), est.P(:,:,i),  est.P(:,:,i+1), est.PP(:,:,i+1), F, Q);
   est.errorsm(:,i)=sqrt((est.xsmooth(:,i)-x(:,i)).^2);
end

    function H=ekfhfun(x,station)
        x1=x(1);y1=x(2);
        x2=station(1);y2=station(2);
        H=[(2*x1 - 2*x2)/(2*((x1 - x2)^2 + (y1 - y2)^2)^(1/2)),0,...
            (2*y1 - 2*y2)/(2*((x1 - x2)^2 + (y1 - y2)^2)^(1/2)),0;
            -(y1 - y2)/((x1 - x2)^2*((y1 - y2)^2/(x1 - x2)^2 + 1)),0,...
            1/((x1 - x2)*((y1 - y2)^2/(x1 - x2)^2 + 1)),0];
    end
end