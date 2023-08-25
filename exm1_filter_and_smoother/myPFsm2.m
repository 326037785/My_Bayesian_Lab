function est=myPFsm2(model,truth,meas)
%% 初值
node=model.xdim;
k=model.K;
F=model.F2;
x=truth.X;R=diag([30,pi/180]).^2;
% a=max(x(1,1:4))+10;b=max(x(3,1:4))+10;
est.X=zeros(model.xdim,k);
est.error=zeros(model.xdim,k);
Q=model.Q;
% xstart=[randi([round(a-a/16),round(a+a/16)],1),20,randi([round(b-b/16),round(b+b/16)],1),20]';
xstart=[meas.Z(1,1),0,meas.Z(2,1),0]';
P=diag([30,3,30,3])^2;
num=model.N;
W=zeros(1,num);
zpre=zeros(model.zdim,num);
xpf=repmat(xstart,[1,num])+10*sqrt(P)*randn(size(xstart,1),num);%随机采样生成粒子
estimatex=zeros(node,1);
Wsm=zeros(num,k);
prob_trans=zeros(num,num);
Wf=zeros(num,k);Wsm=zeros(num,k);xpfilt=zeros(node,num,k);%xsm=zeros(node,num,k);
for i=1:k
    %% 从转移密度中采样
    for j=1:num
        xpf(:,j)=F*xpf(:,j)+1*mvnrnd(zeros(model.xdim,1),Q)';
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
    if ~mod(i,5)
        offset=5;
        if i<10
            offset=4;
        end
        Wsm(:,i)=Wf(:,i);
        for i1=i-1:-1:i-offset
        xtemp=F*xpfilt(:,:,i1)+mvnrnd(zeros(model.xdim,1),Q)';
     for j=1:num
        prob_trans(:,j)=mvnpdf(xpfilt(:,:,i1+1)',xtemp',Q)*Wf(j,i1);
     end
         prob_trans(prob_trans==0)=1e-12;
    prob_trans=prob_trans./(sum(prob_trans)+eps);
%后向平滑
   prob_trans=prob_trans.* Wsm(:,i1+1);
   Wsm(:,i1)=sum(prob_trans,1)';
   % 重采样
%     C= cumsum(Wsm(:,i1));%对权重按列累加求和；
%     index = zeros(1, num); %生成一个100个0的空矩阵；
%     U = linspace(0,1-1/num,num);%生成一个均匀分布U；
%     l=1;
%     for j1=1:num
%         while U(j1)>C(l)
%             l=l+1;
%         end
%         index(j1)=l;
%     end
est.X(:,i1)=xpfilt(:,:,i1)*Wsm(:,i1);
est.error(:,i-1)=sqrt((est.X(:,i1)-x(:,i1)).^2);
        end
    end
end
