function [Uhat, Bhat] = AltGDmin(Tmax,U, Y)

global n1 n2  q r mk n
n=n1*n2;
Uhat = U;
for t=1:1:Tmax
    AU=Afft(Uhat);
    for k =  1 : q
        Bhat(:,k)=(AU(1:mk(k),:,k) ) \ Y(1:mk(k), k);
    end
    Xhat = Uhat * Bhat;
    Grad_U = zeros(n, r);
    
    x_n=zeros(n,q);
    for k=1:1:q
        Yhat(:,k)=AU(:,:,k)*Bhat(:,k);
    end
    Inter=Yhat-Y;
    tmp3=Att(Inter);
    
    Grad_U=tmp3(:,1:q)*(Bhat');
    %%%%%%%%%%%%%%%%% Define eta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if t==1
        w=norm(Grad_U);
        eta=1/(7*w);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Uhat_t0=Uhat;
    Uhat = Uhat - eta * Grad_U;
    [Qu,~]  =  qr(Uhat,0);
    Uhat  =  Qu(:, 1:r);
    %%%%%%%%%%%%%%%%% Loop Exist Condition %%%%%%%%%%%%%%%%%%%%%%%%%%
    Uhat_t1=Uhat;
    Subspace_d= ( norm((Uhat_t0 - Uhat_t1*(Uhat_t1'*Uhat_t0)), 'fro')/sqrt(r));
    if  (Subspace_d <=.01)
        break;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
Yvec_all=Afft(Uhat);
for k =  1 : q
    Bhat(:, k) = (Yvec_all (1:mk(k),:,k)) \ Y(1:mk(k), k);
end
end