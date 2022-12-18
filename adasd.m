clear all
clf
set(0,'defaultaxesfontsize',22);
M=12;
Nb=32;
P=3;
%angles=[25 80 135]*(pi/180);
angdeg=[80 90 115];
angles=angdeg*(pi/180);
dlambda=0.5;
%generate random bits of information
poles=zeros(1,P);
X=zeros(M,Nb);
Rideal=zeros(M,M);
for k=1:P,
mu=pi*cos(angles(1,k));
poles(1,k)=exp(j*mu);
a=exp(j*mu*(0:M-1)).';
Rideal=Rideal+a*a';
br=ones(1,Nb);
temp=rand(1,Nb);
br(find(temp<.5))=-1;
bi=ones(1,Nb);
temp=rand(1,Nb);
bi(find(temp<.5))=-1;
b=br+j*bi;
X=X+a*b;
end
%add some noise
X=X+0.6*(randn(M,Nb)+j*randn(M,Nb));
Rxx=X*X'/Nb;
[E,D,V]=svd(Rxx);
%ESPRIT algorithm:
Es=E(:,1:P);
Es1=Es(1:M-1,:); Es2=Es(2:M,:);
Psi=Es1\Es2;
[T,Phi]=eig(Psi);
Phivec=diag(Phi);
%plot eigenvalues from ESPRIT and compare with true frequencies
polar(0,1,'.')
hold on
plot(real(poles),imag(poles),'kx','MarkerSize',12,'Linewidth',2);
plot(real(Phivec),imag(Phivec),'ro','MarkerSize',12,'Linewidth',2);
hold off
legend('True "poles"','ESPRIT eigenvalues')
