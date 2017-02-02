function xprime = funcD_RV_lag(t,x,z);
global n frprime fvprime rNoise mNoise nScale

% states and lagged states
r = x(1:n); v = x(n+1);
rlag = z(1:n); vlag = z(n+1);

% noise vectors
nR = rNoise(:,ceil((t+1)*nScale));
nM = mNoise(:,ceil((t+1)*nScale));

% derivatives
rprime = frprime(r+nR,v*ones(n,1)+nM);
vprime = fvprime(v,rlag);
xprime = [rprime; vprime];

