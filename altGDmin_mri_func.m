function [Xhat_mri,Time_mri]=altGDmin_mri_func(kdata,b1c,samp)
global nx ny nt nc  b1c samp r mk

 tic;

 y=kdata;
[nx,ny,nt,nc]=size(kdata);

%%%%%%%%% Extracting Mean Image %%%%%%%%%%%%%%%
[zbar_hat,flag,resNE,iter] = cgls_mean(@E_forw_for_mean,@E_back_for_mean, y,0,1e-36,10);
ybar_hat=E_forw_for_mean(zbar_hat);
yinter=y-ybar_hat;
%%%%%%%% Calling altGDMin : Low rank component %%%%%%%%%%%%%%%%%%%
X=altGDmin_basic(yinter);
yk=y-(E_forw(X+zbar_hat));
%%%%%%%%%% Calculating Sparse component %%%%%%%%%%%%
global k
for k=1:nt
    Ehat(:,:,k)=cgls_modi(@E_forw_for_Ak,@E_back_for_Ak,squeeze(yk(:,:,k,:)),0,1e-36,3);
end
Xhat_mri=X+zbar_hat+Ehat;
Time_mri=toc;
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


function Ak = E_forw_for_Ak(xk)
global nx ny nt nc r b1c samp k
s = zeros(nx,ny,nc); 

% smaps_k = squeeze(arg.smaps(:,:,k,:)); %nx ny nc 
samp_k = squeeze(samp(:,:,k)) ;  %nx ny 
s = bsxfun(@times,xk,b1c); % nx ny  nc

S=fft2c_mri(s); % nx ny 1 nc  
Ak = bsxfun(@times,S,samp_k);
end


function x = E_back_for_Ak(zk)
global nx ny nt nc r b1c samp k
zk_2 = squeeze(zk); %zk: nx ny nc
x = zeros(nx,ny); 
% smaps_k = squeeze(arg.smaps(:,:,k,:)); %nx ny nc 
samp_k = samp(:,:,k) ;  %nx ny 

 ztmp = bsxfun(@times,zk_2,samp_k); % nx ny nc 
s = ifft2c_mri(ztmp); %nx ny nc 
x = sum(bsxfun(@times,s,conj(b1c)),3);
end

