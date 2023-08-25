function z=hfun(stationxy,x)
p=size(x,2);
for i=1:p
    z(1,i)=sqrt((stationxy(1)-x(1,i))^2+(stationxy(2)-x(2,i))^2);
    z(2,i)=atan2(x(2,i)-(stationxy(2)),(x(1,i)-stationxy(1)));%弧度
end
