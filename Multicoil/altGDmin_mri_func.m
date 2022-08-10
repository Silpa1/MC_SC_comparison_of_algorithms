function [Xhat_mri,Time_mri]=altGDmin_mri_func(kdata,b1c,samp)
global nx ny nt nc  b1c samp r mk

 tic;

 y=kdata;
[nx,ny,nt,nc]=size(kdata);
param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);

% mask=kdata(:,:,:,1)~=0;
% mask1=repmat(mask,[1,1,1,nc]);


 n=nx*ny;
[zbar_hat,flag,resNE,iter] = cgls_mean(@E_forw_for_mean,@E_back_for_mean, y,0,1e-36,10);
ybar_hat=E_forw_for_mean(zbar_hat);

yinter=y-ybar_hat;
sum_mk=nnz(y);
for k=1:1:nt
    mk(k)=nnz(y(:,:,k,:));
end

m=max(mk);

C_tilda=36;
alpha=C_tilda*norm(yinter(:))^2/sum_mk;
Y_trunc=yinter;
Y_trunc(abs(yinter)>sqrt(alpha))=0;
X0_temp=param.E'*Y_trunc;
DiagMat=diag(ones(1,nt)./sqrt(mk*m));
X0_image=X0_temp;
X0=(reshape(X0_image,[nx*ny,nt]))*(DiagMat);
r_big=floor(min([n/10,nt/10,m/10]));
[Utemp, Stemp,~]=svds(X0,r_big);
SS=diag(Stemp);
E=sum(SS.^2);
Esum=0;
for i=1:1:r_big
    Esum=Esum+((SS(i))^2);
    if Esum >(E*0.85)
        break
    end
end

r=i+1;
r=min(r,r_big);
U0=Utemp(:,1:r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Init End %%%%%%%%%%%%%%%%%%%%%%%
Uhat=U0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% altGDmin Start %%%%%%%%%%%%%%%%%%
T=70;
y_temp=reshape(yinter,[nx*ny,nt,nc]);

for t = 1 : T
    Uhatm=reshape(Uhat,[nx,ny,r]);
    B = E_forw_for_AU_new(Uhatm,y_temp);
    X=reshape(Uhat*B,[nx,ny,nt]);
    Z=param.E'*((param.E*X)-yinter);
    Z_mat=reshape(Z,[nx*ny,nt]);
    Grad_U=Z_mat*B';
    if t==1
        eta=1/(7*norm(Grad_U));
    end
    Uhat_t0=Uhat;
    Uhat=Uhat-eta*Grad_U;
    [Qu,~]  =  qr(Uhat,0);
    Uhat  =  Qu(:, 1:r);
    Uhat_t1=Uhat;
    Subspace_d= ( norm((Uhat_t0 - Uhat_t1*(Uhat_t1'*Uhat_t0)), 'fro')/sqrt(r));
    if  (Subspace_d <=.01)
        break;
    end
end
yk=y-(param.E*(X+zbar_hat));
global k
for k=1:nt
    Ehat(:,:,k)=cgls_modi(@E_forw_for_Ak,@E_back_for_Ak,squeeze(yk(:,:,k,:)),0,1e-36,3);
end
Xhat_mri=X+zbar_hat+Ehat;
Time_mri=toc;
end


function B = E_forw_for_AU_new(U_im,y_temp)
global nx ny nt nc r b1c samp mk
B=zeros(r,nt);
n=nx*ny;
smaps = b1c;
mask1=samp;
% U_im = reshape(Uhat,[nx,ny,r]);
U_im2 = reshape(U_im,[ nx, ny, 1,r]);  % nx ny 1 r  
smaps = reshape(smaps, [nx,ny,nc,1]); % nx ny nc 1
AUtmp = bsxfun(@times, U_im2, smaps);  % nx ny nc r
AUtmp1= fft2c_mri(AUtmp); %nx ny nc r
AUtmp2=reshape(AUtmp1,[nx*ny,nc,r]);
    for k=1:1:nt
        ObsFreq=find(mask1(:,:,k));
        m1 = mk(k);
        AUk_tmp = AUtmp2(ObsFreq,:,: ) ;
        yk_tmp=y_temp(ObsFreq,k,: ) ;
        AUk = reshape(AUk_tmp, [m1, r]);
        yk = reshape(yk_tmp, [m1,1]);
        B(:,k)=AUk\yk;
    end
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
