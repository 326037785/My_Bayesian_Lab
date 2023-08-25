function est=myPF(model,truth,meas)
%% 初值
node=model.xdim;
k=model.K;
H=model.H;F=model.F1;
x=truth.X;R=diag([30,pi/180]).^2;
a=max(x(1,1:4))+10;b=max(x(3,1:4))+10;
est.X=zeros(model.xdim,k);
est.error=zeros(model.xdim,k);
Q=model.Q;
xstart=[randi([round(a-a/16),round(a+a/16)],1),20,randi([round(b-b/16),round(b+b/16)],1),20]';
xstart=[meas.Z(1,1),0,meas.Z(2,1),0]';
P=diag([30,2,30,2])^2;
num=model.N;
W=zeros(1,num);
zpre=zeros(model.zdim,num);
xpf=repmat(xstart,[1,num])+10*sqrt(P)*randn(size(xstart,1),num);%随机采样生成粒子
estimatex=zeros(node,1);
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
    W=1/num;
    estimatex=mean(xpf,2);
    est.X(:,i)=estimatex;
    est.error(:,i)=sqrt((estimatex-x(:,i)).^2);
end