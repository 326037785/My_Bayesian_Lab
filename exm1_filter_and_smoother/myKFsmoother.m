function [xsm, Psm, Psm_future]=myKFsmoother(xsm_future, Psm_future, ...
    xfilt, Pfilt,  Pfilt_future, PPfilt_future, F, Q)
xpred = F*xfilt;
Ppred = F*Pfilt*F' + Q; % Ppred = Cov[X(t+1) | t]
J = Pfilt * F' * inv(Ppred); % smoother gain matrix
xsm = xfilt + J*(xsm_future - xpred);
Psm = Pfilt + J*(Psm_future - Ppred)*J';
Psm_future = PPfilt_future + (Psm_future - Pfilt_future)*inv(Pfilt_future)*PPfilt_future;
