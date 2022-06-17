function [U0] = initAltGDMin( Y)
C_tilda = 6;
global  r m n mk q  

Xhat_trunc = zeros(n, q);
sqrt_alpha = ( C_tilda * norm(Y, 'fro') ) / sqrt(m * q);
Y_trunc = Y;
Y_trunc(abs(Y) > sqrt_alpha) = 0;
sum_true=zeros(n,q);
e=eye(q);
tmp=Att(Y_trunc);
for k=1:1:q
    X_0(:,k)=tmp(:,k)/sqrt(mk(k)*m);
end
%%%%%%%%%%% Automated r %%%%%%%%%%%%%%%%%%%
r_big=floor(min([n/10,q/10,m/10]));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U0=Un(:,1:r);

