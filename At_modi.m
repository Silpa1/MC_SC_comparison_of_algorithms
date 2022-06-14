function x = At_modi(y)
global S2 n1 n2 kk mk m n

%m=max(mk);
w = zeros(n,1);
w(S2(1:mk(kk),kk))=y(1:mk(kk));
tmp3 = n*reshape( ifft2( reshape(w, [n1 n2]) ), [n,1]) ;
x= tmp3;

%Y1vec=reshape(X1vec,[n*q,1]);