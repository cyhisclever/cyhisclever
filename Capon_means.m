%M为阵元数 K为快拍数 X为阵列接收到的信号
function [searching_doa,Pcapon]=Capon_means(M,K,X,d,lambda)
searching_doa=-90:0.1:90;%线阵的搜索范围为-90~90度
R=X*X'/K;
iR=inv(R);
for i=1:length(searching_doa)
    a_theta=exp(-1j*(0:M-1)'*2*pi*d*sin(pi*searching_doa(i)/180)/lambda);
    Pcapon(i)=1./abs((a_theta)'*iR*a_theta);
    Pbf(i)=a_theta'*R*a_theta/(a_theta'*a_theta);
end

