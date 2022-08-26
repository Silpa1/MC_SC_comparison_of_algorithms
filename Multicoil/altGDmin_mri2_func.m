function [Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata,b1c,samp)
global nx ny nt nc  b1c samp r mk

 tic;

 y=kdata;
 [nx,ny,nt,nc]=size(kdata);
%%%%%%%%% Extracting Mean Image %%%%%%%%%%%%%%%
[zbar_hat,flag,resNE,iter] = cgls_mean(@E_forw_for_mean,@E_back_for_mean, y,0,1e-36,10);
ybar_hat=E_forw_for_mean(zbar_hat);
yinter=y-ybar_hat;
%%%%%%%% Calling altGDMin %%%%%%%%%%%%%%%%%%%
X=altGDmin_basic(yinter);

%%%%%%%%%% Calculating Sparse component %%%%%%%%%%%%
param.d = y-E_forw(X+zbar_hat);
param.T=getT(nx,ny,nt);
param.nite=10;
param.tol=0.0025;
M=E_back(param.d);
Ehat=zeros(nx,ny,nt);
param.lambda_S=0.001*max(max(abs(M)));
ite=0;
while(1)
    ite=ite+1;
    M0=M;
    Ehat=param.T'*(SoftThresh(param.T*reshape(M,[nx,ny,nt]),param.lambda_S));
    resk=E_forw(reshape(Ehat,[nx,ny,nt]))-param.d;
    M=Ehat-reshape(E_back(resk),[nx,ny,nt]);
%     tmp2=param.T*reshape(Ehat,[nx,ny,nt]);
    if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
end
Xhat_mri2=X+zbar_hat+Ehat;
Time_mri2=toc;
end



function Amat_zbar  = E_forw_for_mean(zbar) %zbar: nx ny
global nx ny nt nc  b1c samp

smaps = b1c; % nx ny nc
samp = reshape(samp, [nx,ny,nt,1]) ; %nx ny nt 1
s = bsxfun(@times, zbar, smaps);  % nx ny nc
S=fft2c_mri(s); %nx ny nc
Spermute = reshape(S, [nx,ny,1,nc]);  % nx ny 1 nc
Amat_zbar = bsxfun(@times,Spermute,samp);  % nx ny nt nc
end


function zbar_hat  = E_back_for_mean(Y_in) %Y_in: nx ny nt nc
global nx ny nt nc  b1c samp
smaps = b1c; % nx ny nc
samp = reshape(samp, [nx,ny,nt,1]) ; %nx ny nt 1

S = bsxfun(@times,Y_in,samp); %nx ny nt nc
s = ifft2c_mri(S);  %nx ny nt nc
zbar_hat = sum(bsxfun(@times,s,reshape(conj(smaps),[nx,ny,1,nc])),4);  %nx ny nt after the sum step
end


function y=SoftThresh(x,p)
y=(abs(x)-p).*sign(x).*(abs(x)>p);
y(isnan(y))=0;
end 


