function [Xhat_GD_mean,Time_GD_mean]=mean_altGDmin_func(kdata,b1c,samp)
global nx ny nt nc  b1c samp

tic;
[nx,ny,nt,nc]=size(kdata);
y=kdata;
%%%%%%%%% Extracting Mean Image %%%%%%%%%%%%%%%
[zbar_hat,flag,resNE,iter] = cgls_mean(@E_forw_for_mean,@E_back_for_mean, y,0,1e-36,10);
ybar_hat=E_forw_for_mean(zbar_hat);
yinter=y-ybar_hat;
%%%%%%%% Calling altGDMin %%%%%%%%%%%%%%%%%%%
X=altGDmin_basic(yinter);

Xhat_GD_mean=X+zbar_hat;
Time_GD_mean=toc;
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

