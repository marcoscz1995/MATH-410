function bates
%%%%price process
%log_S0 is the initial log equity price
%r interest rate
%q dividend rate
%v0 initial volatility
%lam_y_0 equity jump frequency with poisson distribution
%lam_y_1 
%mu_P_y equity jump size mean
%sig_y equity jump size st dev
%gamma_Y martingale measure process
%rho correlation of W1(t) and W2(t)

%%%variance process paramters
%k mean reversion speed of variance
%theta long term variance
%sig_V volatility of variance
%mu_V variance jump size mean (average variance jump size)
%lambda_V0 volatility of variance jumps (variance of variance jumps, its
%constant)

%M number of monte carlo trials
%N number of time steps per path
%dT daily time steps

%log-price process parameters
r = .01;
q = .01;
lambda_Y_0 = 15;
lambda_Y_1 = 40;
mu_Y = -.02;
sig_Y = .02;
log_S0 = log(90);
gamma_Y = 1;
rho= 0;

%variance process parameters
gamma_V = 1;
k = log(100);
theta = .02;
sig_V = .5;
lambda_V_0 = 15;
mu_V = .03;
V0 = .02;

%other variables
M=100;
N=250;
T=30;
dt=T/N;
K=log(100);

%instantaneous variance
%lambda_Y_t = lambda_Y_0 + lambda_Y_1*(vprocess(:,1));

lambda_Y_t = lambda_Y_0 + lambda_Y_1*V0;
%lambda_Q_Y = phi_P_Z_Y_G*lambda_Y_t;

%momement generating functions of the log equity price and variance jump
%size
phi_P_Z_Y_G=exp(mu_Y*gamma_Y + .5*(sig_Y^2)*(gamma_Y^2));
phi_P_Z_V_G=(1-gamma_V*mu_V)^(-1);

%jump convexity correction
phi_P_Z_Y_1 = exp(mu_Y+.5*(sig_Y^2));
epsilon_P = phi_P_Z_Y_1 -1;

%intensitites
lambda_Q_Y = phi_P_Z_Y_G*lambda_Y_t;
lambda_Q_V = phi_P_Z_V_G*lambda_V_0;

mu_Q_Y = mu_Y + gamma_Y*(sig_Y^2);
mu_Q_V = phi_P_Z_V_G*mu_V;

sig_Q_Y = sig_Y;


%processes initialization
xprocess=zeros(M,N); %option price path
vprocess=zeros(M,N);%volatility path
Y_pprocess=zeros(M,N);%poisson process for Y
V_pprocess=zeros(M,N); %possion process for V
J_Y_process=zeros(M,1);%sum of jumps for Y
J_V_process=zeros(M,1); %sum of jumps for V
dW_V=randn(M,1); %GMB process for volatility
dW_Y=rho.*dW_V + sqrt(1-(rho^2)).*randn(M,1);


%initialize the stock and volatility
xprocess(:,1)=log_S0+(r-q-lambda_Q_Y*epsilon_P-0.5*V0)*dt+sqrt(V0).*dW_Y(:,1)+J_Y_process(:,1);
vprocess(:,1)=V0+k*(theta-V0)*dt+sig_V*sqrt(V0*dt).*dW_V+J_V_process(:,1);

%Monte carlo 
for i=2:N
    %GBM
    dW_V=randn(M,1);
    dW_Y=rho.*dW_V + sqrt(1-(rho^2)).*randn(M,1);
    
    %poisson jumps per day for Y and V
    Y_pprocess(:,1)=poissrnd(lambda_Q_Y*dt,M,1);
    V_pprocess(:,1)=poissrnd(lambda_Q_V*dt,M,1);
    
    %number of jumps for Y and V
    for m=1:M
        J_Y_process(m,1)=sum(normrnd(mu_Q_Y,sig_Y,[1,Y_pprocess(m,1)]));
        V_Y_process(m,1)=sum(exprnd(mu_Q_V,[1,V_pprocess(m,1)]));
    end
    vprocess(:,i)=vprocess(:,i-1)+k*(theta-vprocess(:,i-1)).*dt+sig_V*sqrt(abs(vprocess(:,i-1)).*dt).*dW_V+J_V_process(:,1);
    xprocess(:,i)=xprocess(:,i-1)+(r-q-lambda_Q_Y*epsilon_P-0.5.*vprocess(:,i-1)).*dt+sqrt(abs(vprocess(:,i-1))).*dW_Y(:,1)+J_Y_process(:,1);
end
price = exp(-r*T)*mean(max(0,exp(xprocess(:,N))-K))

end



    
    
    
    








