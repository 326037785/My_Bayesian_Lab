function model=genmodel
%% nonlinear model
model.xdim=4;
model.zdim=2;
model.T=1;
model.K=100;
model.turn=pi/180;
model.A=[1 model.T;0 1];
model.F=[model.A,zeros(2);zeros(2),model.A];
model.F1=[1,(sin(model.turn*model.T)/model.turn),0,(-(1-cos(model.turn*model.T))/model.turn);
           0,(cos(model.turn*model.T)),0,(-sin(model.turn*model.T));
      0,((1-cos(model.turn*model.T))/model.turn),1,(sin(model.turn*model.T)/model.turn);     
      0,(sin(model.turn*model.T)),0,(cos(model.turn*model.T))];
model.F2=[1,(sin(-model.turn*model.T)/-model.turn),0,(-(1-cos(-model.turn*model.T))/-model.turn);
      0,cos(-model.turn*model.T),0,(-sin(-model.turn*model.T));
      0,(1-cos(-model.turn*model.T))/-model.turn,1,(sin(-model.turn*model.T)/-model.turn);
      0,sin(-model.turn*model.T),0,cos(-model.turn*model.T)];
model.H=[1 0 0 0;0 0 1 0];
model.Q=diag([2.25,0.01,2.25,1]).^2;
model.R=diag([10,pi/180]).^2;
model.Pd=1;
model.N=800;