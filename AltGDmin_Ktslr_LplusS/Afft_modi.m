function [Y]= Afft_modi(X)
global S2 n1 n2 n3 q r mk m n kk ww

    tmp = reshape( fft2( reshape(X, [n1 n2]) ), [n,1]) ;
    Y=zeros(m,1);
     Y(1:mk(kk))= tmp(S2(1:mk(kk),kk));
     Y=reshape(Y,[m,1]);
    
  

