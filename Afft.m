function [Y]= Afft(X)
global S2 n1 n2 r mk m n q

for k=1:1:length(X(1,:))
    tmp(:,k) = reshape( fft2( reshape(X(:,k), [n1 n2]) ), [n,1]) ;
end
if (length(X(1,:))==q )
    
    Y=zeros(m,q);
    for k=1:1:q
        Y(1:mk(k),k)= tmp(S2(1:mk(k),k),k);
    end
    
elseif (length(X(1,:))==r)
    
    Y=zeros(m,r,q);
    for k=1:1:q
        Y(1:mk(k),:,k)= tmp(S2(1:mk(k),k),:);
    end
else
        Y=zeros(m,q);
        for k=1:1:q
            Y(1:mk(k),k)=  tmp(S2(1:mk(k),k));
        end
        Y=reshape(Y,[m*q,1]); 
end



