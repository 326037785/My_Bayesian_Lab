function meas=genmeas(model,truth)
k=model.K;
R=model.R;
meas.Z=zeros(model.zdim,k);
stationxy=truth.station([1,3],:);
x=truth.X([1,3],:);
for i=1:k
    meas.Z(:,i)=hfun(stationxy(:,i),x(:,i))+sqrtm(R)*randn(2,1);
end