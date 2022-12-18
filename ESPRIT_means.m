function [Pesprit]=ESPRIT_means(M,K,X,d,lambda,source_number,w)
subM=M-1;
X1=X(1:subM,:);%子阵1接受的数据矢量
X2=X(2:(subM+1),:);%子阵2接受的数据矢量
    %对两个子阵的模型进行合并
    X=[X1;X2];
    R=X*X'/K;
    %对R进行奇异值分解
    [U,S,V]=svd(R);
    R=R-S(2*subM,2*subM)*eye(2*subM);
    [U,S,V]=svd(R);
    Us=U(:,1:source_number);
    disp(Us);
    Us1=Us(1:subM,:);
    Us2=Us((subM+1):2*subM,:);

    %按照公式得到旋转不变矩阵M
    E=pinv(Us1)*Us2;
    disp('E');
    disp(E);
    %对得到的旋转不变矩阵进行特征分解
    [V,D]=eig(E);
    disp(D);
    D=(diag(D)).';
    doa=-asin(angle(D)/pi)*180/pi;
    doa=sort(doa);
    Pesprit=doa;