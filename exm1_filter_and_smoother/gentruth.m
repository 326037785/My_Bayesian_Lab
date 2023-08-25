function truth=gentruth(model)
k=model.K;
truth.X=zeros(model.xdim,k);
truth.station=zeros(model.xdim,k);
stationstart=[randi([0,10],1),randi([5,10],1),randi([0,20],1),randi([5,10],1)]';
xstart=[randi([100,200],1),randi([10,30],1),randi([100,200],1),randi([10,30],1)]';
for i=1:k
    truth.X(:,i)=xstart;
    truth.station(:,i)=stationstart;
    stationstart=model.F*stationstart;
    if i<31
    xstart=model.F1*xstart+sqrt(model.Q)*randn(model.xdim,1);
    end
    if i>30 && i<66
        xstart=model.F2*xstart+sqrt(model.Q)*randn(model.xdim,1);
    end

    if i>67
        xstart=model.F1*xstart+sqrt(model.Q)*randn(model.xdim,1);
    end
end