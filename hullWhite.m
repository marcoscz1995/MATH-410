function price = hullWhite(S0,K,T,k,theta,v0,rho,sigma,r,N,M)
%function to calculate the price of a euro call option
%using the Hull-White model for stochasting volatility
%S0: Initial price of the asset
%K: strike
%T: Time to maturity
%k: mean reversion rate
%theta: long run variance
%v0: current variance
%rho:correlation of W1(t) and W2(t)
%sigma:volatility of volatility paramter
%r: interest rate
%N: number of time steps per path
%M: number of paths in the MC simulation

T = T/252; % annualize the T

dt=T/N; % daily time steps

Call=zeros(M,1);%call is a Mx1 vector of zeros ie call will be priced M times

S_m(1)=S0;
V_m(1)=v0;

for i=1:M %M monte carlo trials
 S_m = zeros(N,1);
 V_m = zeros(N,1);
 
 for j=1:N
     dW1=randn;
     dW2=rho*dW1+sqrt(1-rho^2)*randn; %using choleski decomp
     S_m(j+1)=S_m(j)*r*dt+sqrt(V_m(j)*S_m(j))*dWi;
     V_m(j+1)=k*(theta-V_m(j))*dt+sigma*sqrt(V_m(j))*dW2;
 end
 
Call(j)=exp(-r*T)*(max(0,S_m(N)-K));
end

price = mean(C);


     
 
 
 


