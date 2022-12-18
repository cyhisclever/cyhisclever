function [searching_doa,Pmusic]=Music_mean(M,K,X,d,lambda,source_number,w)
searching_doa=-90:0.1:90;%线阵的搜索范围为-90~90度
R=X*X'/K;
iR=inv(R);
[U,S,V]=svd(R);
Un=U(:,source_number+1:M);
Gn=Un*Un';
for i=1:length(searching_doa)
   a_theta=exp(-1j*(0:M-1)'*2*pi*d*sin(pi*searching_doa(i)/180)/lambda);
   Pmusic(i)=a_theta'*a_theta./abs((a_theta)'*Gn*a_theta);
end
RMSE_MUSIC = sqrt( Pmusic/180/2);