function x = Att(y)
global S2 n1 n2 n3 q  mk m n gg ww n4 tt 
ww=1;
m=max(mk);
if ww==1
    pp=q;
else
    pp=n4;
end
% if tt ==1 & ww==2
%    pp=q- n3;
% end
w = zeros(n,pp);
y_mat=reshape(y, [m,pp]);
for k=1:1:pp
    %    if tt ==1 & ww==2
    %        w(S2(1:gg(k),n3+k),k)=y_mat(1:gg(k),k);
    %    else
    if ww==1
        w(S2(1:gg(k),((ww-1)*pp)+k),k)=y_mat(1:gg(k),k);
    else
        w(S2(1:gg(k),(n3+(ww-2)*pp)+k),k)=y_mat(1:gg(k),k);
    end
    tmp3 = n*reshape( ifft2( reshape(w(:,k), [n1 n2]) ), [n,1]) ;
    x(:,k)= tmp3;
end
%Y1vec=reshape(X1vec,[n*q,1]);