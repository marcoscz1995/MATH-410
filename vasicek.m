function vasicek
%function to price call option with stochastic interest rate under the
%Vasicek model
S0 = 1274;
T=50; %time to maturity on option
K=1200;

%initial parameters
N = 100; %steps
M = 1000; %MC paths
dt = T/N;
sigma = .0072;
rho = -.1; %correlation of W1(t) and W2(t)

%interst rate parameters
alpha = .21; %mean conversion speed
gamma = .1; %mean conversion
rSigma = .05; % interest rate sigma
r0 = .05; %initial interst rate


%initializing the vectors
SNow=S0*ones(M,1); % a M x 1 vector, populated by S_0 
SNext=zeros(M,1); % the M x 1 vector of zeros that will be populated with S_t+1

rNow=r0*ones(M,1);
rNext=zeros(M,1);

for i=1:N
    dW1=randn(M,1);
    dW2=rho.*dW1+sqrt(1-rho^2)*randn; % using choleski decomposition
    rNext = (alpha.*(gamma-rNow).*dt+rSigma.*dW2 +rNow);
    SNext=SNow +(dt.*rNow).*SNow+sigma.*SNow.*dW1;
    rNow=rNext;
    SNow=SNext;    
end
 price=mean(exp(-rNow*T)*mean(max(0,SNow-K))) %take the mean because rNow is a vector
end
 

