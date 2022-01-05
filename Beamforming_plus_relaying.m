clc 
clear all
tic
%Parameters of Simulation:

N=3;
Ps=10; %10dB
I1=10^(10/10); %10dB
I2=10^(10/10); %10dB
Pr=0:4:20; %db scale
no_iterations=50;
no_relays=3; %k
Rate_vector=zeros(1,no_relays);
R=zeros(1,no_iterations);
R_ave=zeros(1,length(Pr));

for Pr_index=1:length(Pr)
for iterations=1:no_iterations
for k=1:no_relays


%normizations and Relations:

h_1k=1/sqrt(2)*(randn(N,1)+1i*randn(N,1));
h_2k=1/sqrt(2)*(randn(N,1)+1i*randn(N,1));
g_21=1/sqrt(2)*(randn(N,1)+1i*randn(N,1));
g_22=1/sqrt(2)*(randn(N,1)+1i*randn(N,1));
g_11=1/sqrt(2)*(randn(1,1)+1i*randn(1,1)); %Channel (SU_tx - PU_1)
g_12=1/sqrt(2)*(randn(1,1)+1i*randn(1,1)); %Channel (SU_tx - PU_2)
h=kron(conj(h_1k),h_2k);
H_1=kron(conj(h_1k),eye(N));
H_2=kron(eye(N),h_2k);
G_11=kron(conj(h_1k),g_21);
G_21=kron(eye(N),g_21);
G_12=kron(conj(h_1k),g_22);
G_22=kron(eye(N),g_22);
sigma_2=1;

Power_SU=min(Ps,min((I1/abs(g_11)^2),(I2/abs(g_12)^2)));


A_1=Power_SU*h*(h)';
A_2=sigma_2*H_2*(H_2)';
A_3=Power_SU*H_1*(H_1)'+sigma_2*eye(N^2);
B_1=Power_SU*G_11*(G_11)'+sigma_2*G_21*(G_21)';
B_2=Power_SU*G_12*(G_12)'+sigma_2*G_22*(G_22)';



cvx_begin

variable S(N^2,N^2) hermitian semidefinite;

variable u nonnegative;

maximize real(trace(A_1*S));


subject to 

real(trace(A_2*S)+sigma_2*u) == 1;

real(trace(A_3*S)) <= u* 10^(Pr(Pr_index)/10);

real(trace(B_1*S)) <= u*I1;

real(trace(B_2*S)) <= u*I2;


cvx_end

gamma=cvx_optval;


cvx_begin

variable W(N^2,N^2) hermitian semidefinite;

minimize real(trace(A_3*W));

subject to 

delta_gamma=1e-3;

real(trace(A_1*W))-(gamma-delta_gamma) * real(trace(A_2*W)+sigma_2)>=0;

real(trace(B_1*W)) <= I1;

real(trace(B_2*W)) <= I2;

cvx_end

W;

[V,D]=eig(W);

f=sqrt(D(end,end))*V(:,end);

F=reshape(f,N,N);

Rate_vector(k)=0.5*log2(1+(Power_SU*(abs(h_2k'*F*h_1k))^2)...
    /((sigma_2*(norm(h_2k'*F))^2+sigma_2)));

end

R(iterations)= max(Rate_vector);

end

R_ave(Pr_index)=sum(R)/iterations;

end

%Plot:

plot(Pr,R_ave)
xlabel ('P_{R}(dB)');
ylabel('R_{ave}(bps/Hz)');

toc
