function [Uhat] = initAltGDMin( Y)
Cy = 6;
%[m,q] = size(S);
global S n1 n2  r m n mk ww q gg r_big

[nn,n3]=size(Y);
Xhat_trunc = zeros(n, n3);
threshold = ( Cy * norm(Y, 'fro') ) / sqrt(m * n3);
Y_trunc = Y;
Y_trunc(abs(Y) > threshold) = 0;
sum_true=zeros(n,n3);
e=eye(n3);
tmp=Att(Y_trunc);
% tic;
% for k=1:1:n3
%     sum_true=sum_true+(tmp(:,k)*(e(:,k)')/gg(k));
% end
% t1=toc;
% tic;
for k=1:1:n3
    X_0(:,k)=tmp(:,k)/gg(k);
end
% t2=toc;
AA=[n/10,n3/10];
    r_big=floor(min(AA));
    
  % [Xbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,20);
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
%             r=min(r,r_big);
if r>r_big
    [Un,Sigman,Vn]=svds(double(X_0),r);
end
    
Uhat=Un(:,1:r);

%%%%%%%%%%%%%%%%%%%%%Finding Rank%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=-diag(Shat);
% diffs= diff(A);
% result=max( diffs(diffs>=0) );
% find(diffs==result);
% rank_1=find(diffs==result);
% S_tmp=diag(Shat);
% K=sum(S_tmp);
% B=0.9*K;
% P=S_tmp(1);
% for i=2:1:length(S_tmp)
%     P=P+S_tmp(i);
%     if (P >=B)
%         break;
%     end
% end
% r_hat=i;
% A=[n/10,q/10,r_hat];
% r=min(A);
% Uhat=Uhat(:,1:r);
% for i=1:1:(length(S_tmp)-1)
%     L(i)=diffs(i)/S_tmp(i);    
% end
% L=L';
% result1=max( L(L>=0) )
% rank_3=find(L==result1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%One dimensional code
% % for k = 1 : q
% %     Xhat_trunc(:, k) = sum( A(:,:,k)' * diag(Y_trunc(:, k)) , 2 );
% % end
% % [Uhat1, Shat1, ~] = svds(Xhat_trunc , r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%One dimensional code End
rrrr=2;