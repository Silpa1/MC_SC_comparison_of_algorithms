function [Uhat] = initAltGDMin( Y)
Cy = 6;
%[m,q] = size(S);
global  r m n mk q  

Xhat_trunc = zeros(n, q);
threshold = ( Cy * norm(Y, 'fro') ) / sqrt(m * q);
Y_trunc = Y;
Y_trunc(abs(Y) > threshold) = 0;
sum_true=zeros(n,q);
e=eye(q);
tmp=Att(Y_trunc);
% tic;
% for k=1:1:n3
%     sum_true=sum_true+(tmp(:,k)*(e(:,k)')/gg(k));
% end
% t1=toc;
% tic;
for k=1:1:q
    X_0(:,k)=tmp(:,k)/mk(k);
end
% t2=toc;
AA=[n/10,q/10,m/10];
r_big=floor(min(AA));

[Un,Sigman,Vn]=svds(double(X_0),r_big);
SS=diag(Sigman);
E=sum(SS.^2);
Esum=0;
for i=1:1:r_big
    Esum=Esum+((SS(i))^2);
    if Esum >(E*0.85)
        break
    end
end
r=i+1;
r=min(r,r_big);
Uhat=Un(:,1:r);

