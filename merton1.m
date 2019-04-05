function jumps
% this simulates a stochastic process with jumps
%S0 initital stock pirce
%r interest rate
%d dividend rate
%sigma asset volatility
%lambda is a historical value of the yearly average number of jumps
%muJ: jump mean parameter
%sigmaJ: jump volatiity
%N: number of time steps per year
%M: number of monte carlo trials

N=252;
M=1000;
muJ=0;
sigmaJ=.3411;
sigma=.3;
r=.01;
S0=1274;
lambda=3;
T=50;
K=1200;
T = T/252; % annualize the T

dt=T/N; % daily time steps

%initialize the vectors

SNow=S0*ones(M,1); % a M x 1 vector, populated by S_0 
SNext=zeros(M,1); % the M x 1 vector of zeros that will be populated with S_t+1


for i=1:N
    dN=poissrnd(lambda*dt,M,1); %number of jumps in a given day
    dW=randn(M,1); %GMB
    dZ=randn(M,1);
    logY=(muJ)*exp(-.5*sigmaJ+sqrt(sigmaJ).*dZ); %jump size distribution
    SNext=SNow+(dt*r-lambda*dt).*SNow+sigma.*SNow.*dW+SNow.*(logY-1).*dN; %need to correct for bias so we remove the lambda*dt from SNow
    SNow=SNext;
end
price=exp(-r*T)*mean(max(0,SNow-K))
end



    





