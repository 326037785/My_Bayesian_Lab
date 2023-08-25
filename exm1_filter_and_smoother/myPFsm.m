function est=myPFsm(model,truth,meas)
%% 初值
node=model.xdim;
k=model.K;
H=model.H;F=model.F2;
x=truth.X;R=diag([30,pi/180]).^2;
a=max(x(1,1:4))+10;b=max(x(3,1:4))+10;
est.X=zeros(model.xdim,k);
est.error=zeros(model.xdim,k);
Q=model.Q;
xstart=[randi([round(a-a/16),round(a+a/16)],1),20,randi([round(b-b/16),round(b+b/16)],1),20]';
xstart=[meas.Z(1,1),0,meas.Z(2,1),0]';
P=diag([30,3,30,3])^2;
num=model.N;
W=zeros(1,num);
zpre=zeros(model.zdim,num);
xpf=repmat(xstart,[1,num])+sqrt(P)*randn(size(xstart,1),num);%随机采样生成粒子
estimatex=zeros(node,1);
Wsm=W;
Wf=zeros(num,k);Wsm=zeros(num,k);xpfilt=zeros(node,num,k);xsm=zeros(node,num,k);
for i=1:k
    %% 从转移密度中采样
    for j=1:num
        xpf(:,j)=F*xpf(:,j)+2*mvnrnd(zeros(model.xdim,1),Q)';
        zpre(:,j)=hfun(truth.station([1,3],i),xpf([1,3],j));
        W(j) = mvnpdf(meas.Z(:,i),zpre(:,j),R)+eps;%计算权重
%         W(j) = inv(sqrt(2*pi*det(R)))*exp(-.5*(zpre(:,j)-meas.Z(:,i))'*inv(R)...
%         *(zpre(:,j)-meas.Z(:,i)))+eps;
    end
    W = W./(sum(W));
    C= cumsum(W);%对权重按列累加求和；
    index = zeros(1, num); %生成一个100个0的空矩阵；
    U = linspace(0,1-1/num,num);%生成一个均匀分布U；
    l=1;
    for j=1:num
        while U(j)>C(l)
            l=l+1;
        end
        index(j)=l;
    end
    xpf=xpf(:,index);
    xpfilt(:,:,i)=xpf;
    W=1/num;
    Wf(:,i)=W(:);%或者取重采样之后
    estimatex=mean(xpf,2);
    est.X(:,i)=estimatex;
    est.error(:,i)=sqrt((estimatex-x(:,i)).^2);
end
%% 后向平滑
est.XSM=zeros(size(est.X));
est.errorsm=est.error;
est.XSM(:,end)=est.X(:,end);
Wsm(:,end)=Wf(:,end);
prob_trans=zeros(num,num);
for i=k-1:-1:1
    for j=1:num
        xtemp=F*xpfilt(:,j,i)+2*mvnrnd(zeros(model.xdim,1),Q)';
        prob_trans(:,j)=mvnpdf(xpfilt(:,:,i+1)',xtemp',Q)*Wf(j,i);
    end
    prob_trans(prob_trans==0)=1e-12;
    prob_trans=prob_trans./(sum(prob_trans)+eps);
%后向平滑
   prob_trans=prob_trans.* Wsm(:,i+1);
   Wsm(:,i)=sum(prob_trans,1)';
   % 重采样
    C= cumsum(Wsm(:,i));%对权重按列累加求和；
    index = zeros(1, num); %生成一个100个0的空矩阵；
    U = linspace(0,1-1/num,num);%生成一个均匀分布U；
    l=1;
    for j1=1:num
        while U(j1)>C(l)
            l=l+1;
        end
        index(j1)=l;
    end
    %状态提取
   est.XSM(:,i)=xpfilt(:,:,i)*Wsm(:,i);
   est.errorsm(:,i)=sqrt((est.XSM(:,i)-x(:,i)).^2);
end