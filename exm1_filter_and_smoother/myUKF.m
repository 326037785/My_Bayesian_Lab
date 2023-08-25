    function est=myUKF(model,truth,meas,opt)
%% 初值赋予
node=model.xdim;
k=model.K;
F=model.F1;
x=truth.X;R=model.R;
a=max(x(1,1:4))+10;b=max(x(3,1:4))+10;
est.X=zeros(model.xdim,k);
est.error=zeros(model.xdim,k);
Q=model.Q;
xstart=[randi([round(a-a/16),round(a+a/16)],1),20,randi([round(b-b/16),round(b+b/16)],1),20]';
P=diag([10,2,10,2])^2;
Wc=zeros(1,2*node+1);
Wm=zeros(1,2*node+1);
zsig=zeros(model.zdim,2*node+1);
xsig=zeros(model.xdim,2*node+1);
P_ukf=zeros(model.xdim,model.xdim,k);
switch opt
    case 1
        %% UKF参数分配
        beta=2;%用于描述状态的分布情况，高斯最佳为2
        alpha=0.01;%小的正值，描述x的分散程度
        tao=0;
        lamda=alpha*(node+tao)-1;
        %% 均值和协方差权重
        for i=2:2*node+1
            Wm(i)=0.5/(node+lamda);
            Wc(i)=0.5/(node+lamda);
        end
        Wm(1)=lamda/(node+lamda);
        Wc(1)=lamda/(node+lamda)+(1-alpha^2+beta);
        %% 滤波
        for i=1:k
            cho=chol((node+lamda)*P);
            xsig(:,1)=xstart;
            for j=1:node
                xsig(:,j+1)=xstart+cho(:,j);
                xsig(:,node+1+j)=xstart-cho(:,j);
            end
            xsignew=F*xsig;
            xstart=sum(xsignew.*Wm,2);
            tmp=repmat(xstart,[1,2*node+1]);
            P=(Wc.*(tmp-xsig)*(tmp-xsig)')+Q;
            P=(P+P')/2;
            cho2=chol((node+lamda)*P)';
            xsig(:,1)=xstart;
            for j=1:node
                xsig(:,j+1)=xstart+cho2(:,j);
                xsig(:,node+1+j)=xstart-cho2(:,j);
            end
            %% 求预测的sigma量测
            zsig=hfun(truth.station([1,3],i),xsig([1,3],:));
            zpre=sum(zsig.*Wm,2);
            tmp2=repmat(zpre,[1,2*node+1]);
            Pzz=(Wc.*(tmp2-zsig)*(tmp2-zsig)')+R;
            Pxz=(Wc.*(tmp-xsig)*(tmp2-zsig)');
            K=Pxz*inv(Pzz);
            xstart=xstart+K*(meas.Z(:,i)-zpre);
            P=P-K*Pzz*K';
            P_ukf(:,:,i)=P;
            est.X(:,i)=xstart;
            est.error(:,i)=sqrt((xstart-x(:,i)).^2);
        end
        %% 平滑
        flag=1;
        if flag==1
            xsig=zeros(model.xdim,2*node+1);
            est.Xsm=zeros(model.xdim,k);est.Psm=zeros(model.xdim,model.xdim,k);
            est.Xsm(:,end)=est.X(:,end);est.Psm(:,:,end)=P_ukf(:,:,end);
            est.errorsm=zeros(model.xdim,k);est.errorsm(:,end)=est.error(:,end);
            for i=k-1:-1:1
                cho=chol((node+lamda)*P_ukf(:,:,i));
                xstart=est.X(:,i);
                xsig(:,1)=xstart;
                xsig(:,2:node+1)=xstart+cho(:,1:node);
                xsig(:,node+2:end)=xstart-cho(:,1:node);
                xsignew=F*xsig;
                xstart=sum(xsignew.*Wm,2);
                tmp=repmat(xstart,[1,2*node+1]);
                tmp2=repmat(est.X(:,i),[1,2*node+1]);
                P=(Wc.*(tmp-xsignew)*(tmp-xsignew)');
%                 cho2=chol((node+lamda)*P);
%                 xsignew(:,2:node+1)=xstart+cho2(:,1:node);
%                 xsignew(:,node+2:end)=xstart-cho2(:,1:node);
                P=(P+P')/2;
                C=(Wc.*(tmp2-xsig)*(tmp-xsignew)');
                D=C/(P);
                xstart=est.X(:,i)+D*(est.Xsm(:,i+1)-xstart);
                est.Xsm(:,i)=xstart;
                est.Psm(:,:,i)=P_ukf(:,:,i)+D*(est.Psm(:,:,i+1)-P)*D';
                est.errorsm(:,i)=sqrt((xstart-x(:,i)).^2);
            end

        end




    case 2
        %% UKF参数分配
        beta=2;%用于描述状态的分布情况，高斯最佳为2
        alpha=1;%小的正值，描述x的分散程度
        lamda=3-node;
        %% 均值和协方差权重
        for i=2:2*node+1
            Wm(i)=0.5/(node+lamda);
            Wc(i)=0.5/(node+lamda);
        end
        Wm(1)=lamda/(node+lamda);
        Wc(1)=lamda/((node+lamda)+(1-alpha^2+beta));
        %% 滤波
        for i=1:k
            cho=chol((node+lamda)*P);
            xsig(:,1)=xstart;
            for j=1:node
                xsig(:,j+1)=xstart+cho(:,j);
                xsig(:,node+1+j)=xstart-cho(:,j);
            end
            xsignew=F*xsig;
            xstart=sum(xsignew.*Wm,2);
            tmp=repmat(xstart,[1,2*node+1]);
            P=(Wc.*(tmp-xsig)*(tmp-xsig)')+Q;
            P=(P+P')/2;
            cho2=chol((node+lamda)*P)';
            xsig(:,1)=xstart;
            for j=1:node
                xsig(:,j+1)=xstart+cho2(:,j);
                xsig(:,node+1+j)=xstart-cho2(:,j);
            end
            %% 求预测的sigma量测
            zsig=hfun(truth.station([1,3],i),xsig([1,3],:));
            zpre=sum(zsig.*Wm,2);
            tmp2=repmat(zpre,[1,2*node+1]);
            Pzz=(Wc.*(tmp2-zsig)*(tmp2-zsig)')+R;
            Pxz=(Wc.*(tmp-xsig)*(tmp2-zsig)');
            K=Pxz*inv(Pzz);
            xstart=xstart+K*(meas.Z(:,i)-zpre);
            P=P-K*Pzz*K';
            P_ukf(:,:,i)=P;
            est.X(:,i)=xstart;
            est.error(:,i)=sqrt((xstart-x(:,i)).^2);
        end
        %% 平滑
        flag=1;
        if flag==1
            xsig=zeros(model.xdim,2*node+1);
            est.Xsm=zeros(model.xdim,k);est.Psm=zeros(model.xdim,model.xdim,k);
            est.Xsm(:,end)=est.X(:,end);est.Psm(:,:,end)=P_ukf(:,:,end);
            est.errorsm=zeros(model.xdim,k);est.errorsm(:,end)=est.error(:,end);
            for i=k-1:-1:1
                cho=chol((node+lamda)*P_ukf(:,:,i));
                xstart=est.X(:,i);
                xsig(:,1)=xstart;
                xsig(:,2:node+1)=xstart+cho(:,1:node);
                xsig(:,node+2:end)=xstart-cho(:,1:node);
                xsignew=F*xsig;
                xstart=sum(xsignew.*Wm,2);
                tmp=repmat(xstart,[1,2*node+1]);
                tmp2=repmat(est.X(:,k),[1,2*node+1]);
                P=(Wc.*(tmp2-xsignew)*(tmp-xsignew)')+Q;
%                 cho2=chol((node+lamda)*P);
%                 xsignew(:,2:node+1)=xstart+cho2(:,1:node);
%                 xsignew(:,node+2:end)=xstart-cho2(:,1:node);
%                 P=(P+P')/2;
                C=P_ukf(:,:,i+1)*F*inv(P);
%                 C=(Wc.*(tmp2-xsig)*(tmp-xsignew)');
                D=C/P;
                xstart=est.X(:,i)+D*(est.Xsm(:,i+1)-xstart);
                est.Xsm(:,i)=xstart;
                est.Psm(:,:,i)=P_ukf(:,:,i)+D*(est.Psm(:,:,i+1)-P)*D';
                est.errorsm(:,i)=sqrt((xstart-x(:,i)).^2);
            end

        end

end

