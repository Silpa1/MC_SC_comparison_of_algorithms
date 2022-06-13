function [Uhat, Bhat] = GDMin_wi(T,U, Y)

global S2 n1 n2  q r mk m n gg n3 
[nn,qq]=size(Y);
ErrorU = [];
ErrorUF = [];
ErrorX = [];
ExeTime = [];
%tic;
%[U_init] = initAltGDMin(Y);
n=n1*n2;
%[m, q] = size(S);
Uhat = U;
%[Usta, Sstar, Vstar]=svd (X);
%Ustar= Usta(:,1:r);
% ErrorU = [abs(sin(subspace(Ustar, Uhat))), ErrorU];
%ErrorUF(1) =( norm((Uhat - Ustar*(Ustar'*Uhat)), 'fro')/sqrt(r));
%ExeTime = [toc, ExeTime];
uu=1;
for t=1:1:T
    % Update of B
    %tmp3=A_Cardiac_U(Uhat);
    
   AU=Afft(Uhat);
    for k =  1 : n3
        
        %Bhat(:, k) = (tmp3(:,:,k) ) \ Y(:, k);
        Bhat(:,k)=(AU(1:gg(k),:,k) ) \ Y(1:gg(k), k);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%One dimensional code
        %% Bhat1(:, k) = ( A(:, :, k) * Uhat ) \ Y(:, k);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%One dimensional code End
    end
    %norm(Bhat1-Bhat)
    Xhat = Uhat * Bhat;
    % Update U
    Grad_U = zeros(n, r);

     x_n=zeros(n,n3);
     for k=1:1:n3
        Yhat(:,k)=AU(:,:,k)*Bhat(:,k);
        %Inter(1:gg(k), k) = Yhat(1:gg(k),k)- Y(1:gg(k), k);
     end
    Inter=Yhat-Y;
    tmp3=Att(Inter);
%     tic;
%     for k =1:1:qq
%         Grad_U = Grad_U+(tmp3(:,k)*(Bhat(:,k)'));
%     end
%     t1=toc;
    Grad_U=tmp3(:,1:n3)*(Bhat');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%One dimensional code
        %% Grad1_U = Grad1_U + A(:, :, k)' * (Inter(:, k)) * Bhat(:, k)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%One dimensional code End
        if t==1
            w=norm(Grad_U);
            eta=1/(7*w);
        end
        Uhat_t0=Uhat;
        Uhat = Uhat - eta * Grad_U;
        [Qu,~]  =  qr(Uhat,0);
        Uhat  =  Qu(:, 1:r);
        Uhat_t1=Uhat;
          Subspace_d= ( norm((Uhat_t0 - Uhat_t1*(Uhat_t1'*Uhat_t0)), 'fro')/sqrt(r));
         if  (Subspace_d <=.01)
             break;
         end 
       uu=  uu+1;

       
    end
    Yvec_all=Afft(Uhat);
    for k =  1 : n3
        Bhat(:, k) = (Yvec_all (1:gg(k),:,k)) \ Y(1:gg(k), k);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%One dimensional code
        %% Bhat1(:, k) = ( A(:, :, k) * Uhat ) \ Y(:, k);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%One dimensional code End
        
    end
end