function [Y]= Afft(X)
global S2 n1 n2 n3 r mk m n ww gg  q
m=max(gg);
% n3=q;
for k=1:1:length(X(1,:))
    tmp(:,k) = reshape( fft2( reshape(X(:,k), [n1 n2]) ), [n,1]) ;
end
if (length(X(1,:))==n3 )
    
    if ww==1
        Y=zeros(m,n3);
        for k=1:1:n3
            Y(1:gg(k),k)= tmp(S2(1:gg(k),k),k);
        end
    else
        Y=zeros(m,n3);
        for k=1:1:n3
            Y(1:gg(k),k)= tmp(S2(1:gg(k),(n3+(ww-2)*n3+k)),k);
        end
    end
elseif (length(X(1,:))==r)
    if ww==1
        Y=zeros(m,r,n3);
        for k=1:1:n3
            Y(1:gg(k),:,k)= tmp(S2(1:gg(k),k),:);
        end
    else
        Y=zeros(m,r,n3);
        for k=1:1:n3
            Y(1:gg(k),:,k)= tmp(S2(1:gg(k),(n3+(ww-2)*n3+k)),:);
        end
        
    end
else
    if ww==1
        Y=zeros(m,n3);
        for k=1:1:n3
            Y(1:gg(k),k)=  tmp(S2(1:gg(k),((ww-1)*n3+k)));
        end
        Y=reshape(Y,[m*n3,1]);
    else
        Y=zeros(m,n3);
        for k=1:1:n3
            Y(1:gg(k),k)=  tmp(S2(1:gg(k),(n3+(ww-2)*n3+k)));
        end
        Y=reshape(Y,[m*n3,1]);
    end
end



