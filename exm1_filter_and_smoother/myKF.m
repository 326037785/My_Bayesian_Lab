function varargout=myKF(varargin)
%default %%默认是4维
%P是协方差，
%例子1：[x,n]=myKF('X',x,'p',P)
% PP是Cov[ X(t), X(t-1) | z(:, 1:t) ]
X=[0;2;0;2];
F=[1 1 0 0;
    0 1 0 0;
    0 0 1 1;
    0 0 0 1];
H=[1 0 0 0;
    0 0 1 0];
P=eye(4)*9;
Q=diag([5;1;5;1]).^2;
R=diag([2;2]).^2;
Z=[9;11];
I=eye(size(X,1));
for k=1:2:length(varargin)
    switch upper(varargin{k})
        case 'X'
            X=varargin{k+1};
        case 'F'
            F=varargin{k+1};
        case 'H'
            H=varargin{k+1};
        case 'P'
            P=varargin{k+1};
        case 'Q'
            Q=varargin{k+1};
        case 'R'
            R=varargin{k+1};
        case 'Z'
            Z=varargin{k+1};
    end
end
X=F*X;
Pre=F*P*F'+Q;
if isempty(Z)
    P=Pre;
    PP=F*P*F';
else
    K=Pre*H'*inv(H*Pre*H'+R);
    X=X+K*(Z-H*X);
    PP=(I-K*H)*F*P;
    P=(I-K*H)*Pre;
    est.X=X;
    est.N=size(X,2);
    est.P=P;
    est.PP=PP;
end
switch nargout
    case 1
    varargout{1}=est;

    case 2
        varargout{1}=est.X;
        varargout{2}=est.N;
    case 3
        varargout{1}=est.X;
        varargout{2}=est.N;
        varargout{3}=est.P;
    case 4
        varargout{1}=est.X;
        varargout{2}=est.N;
        varargout{3}=est.P;
        varargout{4}=est.PP;
    otherwise
    msgbox('你输错了，重来！！！！！！');
end
