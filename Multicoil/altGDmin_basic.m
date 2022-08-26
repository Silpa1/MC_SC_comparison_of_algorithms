function [X]=  altGDmin_basic(y)
global nx ny nt nc n r mk samp b1c
sum_mk=nnz(y);
for k=1:1:nt
    mk(k)=nnz(y(:,:,k,:));
end

m=max(mk);

C_tilda=36;
alpha=C_tilda*norm(y(:))^2/sum_mk;
Y_trunc=y;
Y_trunc(abs(y)>sqrt(alpha))=0;
% X0_temp=param.E'*Y_trunc;
X0_temp=E_back(Y_trunc);
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
y_temp=reshape(y,[nx*ny,nt,nc]);

for t = 1 : T
    Uhatm=reshape(Uhat,[nx,ny,r]);
    B = E_forw_for_AU_new(Uhatm,y_temp);
    X=reshape(Uhat*B,[nx,ny,nt]);
     Z= E_back(E_forw(X)-y);
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