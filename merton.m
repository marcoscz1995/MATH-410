function price = merton(S0,K,T,k,lambda,muJ,sigmaJ,,r,N,M)
%function to calculate the price of a euro call option
%using the Heston model for stochasting volatility
%S0: Initial price of the asset
%K: strike
%T: Time to maturity
%r: interest rate
%sigma: volatility
%k: negative value of the expected percentage jump
%lambda: average number of jumps per year
%muJ: jump mean paramter
%sigmaJ: jump volatility
%N: number of time steps per path
%M: number of paths in the MC simulation

T = T/252; % annualize the T

dt=T/N; % daily time steps

k=exp(-muJ)-1 %negative value of the expected percentage jump

%simulating the poisson random variables
poiss=poissrnd(lambda*dt,1,N);
dZ=randn(1,N);
poissJumps=zeros(1,N);

%determinitation of jumps for the poisson process
for i=1:N
    poissJumps(1,i)=sum(randn(poiss(1,i),1));
end

%constructing the increments
incrementJumps=(r-lambda*k)*dt*ones(1,N)+sigma*sqrt(dt)*dZ+(muJ-.5*sigmaJ


Call=zeros(M,1);%call is a Mx1 vector of zeros ie call will be priced M times

S_m(1)=S0;

for i=1:M %M monte carlo trials
 S_m = zeros(N,1);
 
 for j=1:N
     
     Poiss=poissrnd(lambda*dt,M,1); %poiss dist number of jumps
     
     
     %simulate the number of jumps
     t=0; %time 0
     tau=[]; %vector of the jumps
     while t<T
         dt = -log(rand)/lambdal % jump time
         t = t+dt; %add the jump time
         tau=[tau;dt]; % make a matrix
     end
     tau(end)= T-(t-dt);
     n=length(tau); %number of jumps
     W1 = randn(1,n);
     
     
     
     dW1=randn;
     dW2=rho*dW1+sqrt(1-rho^2)*randn; %using choleski decomp
     S_m(j+1)=S_m(j)*r*dt+sqrt(V_m(j)*S_m(j))*dWi;
     V_m(j+1)=k*(theta-V_m(j))*dt+sigma*sqrt(V_m(j))*dW2;
 end
 
Call(j)=exp(-r*T)*(max(0,S_m(N)-K));
end

price = mean(C);


     
 
 
 


