c=3e8;
fc=1e10;
%lambda=c/fc;
w=[pi/4 pi/6 pi/3].';%信号频率
source_number=length(w);%信元数
lambda=sum(2*pi*3e8./w)/3;%信号波长 
d=0.5*lambda;
alpha=[-17, -6,20];   % True DOAs
SNR=10;
K=1000;                  % The number of snapshots
M=10;                  % The number of sensors
%% 生成信号
N_alpha=length(alpha);
A=exp(-i*pi*(0:M-1)'*sin(alpha*pi/180));
Vj=diag(sqrt((   10.^(SNR/10)   )/2));
S=Vj*(randn(N_alpha,K)+i*randn(N_alpha,K));
noise=sqrt(1/2)*(randn(M,K)+i*randn(M,K));
X=A*S+noise;
%% 稀疏贝叶斯算法
resolution=1;          % grid interval
search_area=[-90:resolution:90];   % grid 
etc=M;                             %   etc may be chosen from 1 to M.
X=signal(M, alpha, SNR, K,w);    
[Pm_root,search_root]=Bayesian_DOA_root(X,search_area,etc);%稀疏贝叶斯法
%% 其余方法
[searching_doa,Pcapon]=Capon_means(M,K,X,d,lambda);%Capon法
[searching_doa,Pbf]=PBF_means(M,K,X,d,lambda);%CBF法
[searching_doa,Pmusic]=Music_mean(M,K,X,d,lambda,source_number,w);%music法
[Pesprit]=ESPRIT_means(M,K,X,d,lambda,source_number,w)%ESPRIT法
%% omp法
phi=i*randn(K,180)+randn(K,180);
[theta] = CS_OMP( Y,phi,20 )
y=ones(1,length(alpha));
%% 画图
figure
plot(searching_doa,10*log10(Pcapon/max(Pcapon)),'k-',search_root,10*log10(Pm_root/max(Pm_root)),'b-',searching_doa,10*log10(Pbf/max(Pbf)),'r-',searching_doa,10*log10(Pmusic/max(Pmusic)),'y-');
grid on;
hold on;
stem(Pesprit,-60*ones(1,length(alpha)),'g-');
legend('Capon Spectrum','Bayes Spectrum','Cbf Spectrum','Music Spectrum','Esprit Sepctrum');
%% 

  