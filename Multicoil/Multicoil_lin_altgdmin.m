function[] = Multicoil_lin_altgdmin(fid, fid2,name,radial,kdata,b1,samp, Xtrue)
global r mk
mk=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ktslr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  [nx,ny,nt,nc]=size(kdata);
% % tmp = sqrt(sum(abs((b1)).^2,3));
% % b1c = div0(b1,tmp);
% % tmp = max(Xtrue(:));
% % Xtrue = div0(Xtrue,tmp);
% % param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);
% % % A = @(z)A_fhp3D_p(z);
% % % At=@(z)At_fhp3D_p(z);
% % b=kdata;
% % tic;
% % x_init = param.E'*kdata;
% % mu =20e-6;
% % opts.mu = mu;
% % opts.betarate = 35;
% % opts.outer_iter =12;
% % opts.inner_iter = 50;
% % opts.beta=1e3; %
% % 
% % U=x_init; [m n d] = size(U);
% % Lam= zeros(m,n,d);
% % 
% % o=0;
% % earray=zeros(1,opts.outer_iter*opts.inner_iter); cost=zeros(1,opts.outer_iter*opts.inner_iter);
% % 
% % fac=length((find(b~=0)));
% % U = double(U); b=double(b); fac=double(fac);
% % for out = single(1:opts.outer_iter),
% %     o=o+1;
% %     for in = single(1:opts.inner_iter)
% %         
% %         z = itemp_FT(U);
% %         
% %         
% %         z1 = (abs(z)-1/opts.beta); z1 = z1.*(z1>0);
% %         z2 = abs(z) + (abs(z)<1e-12);
% %         z = z1.*z./z2;
% %         z = temp_FT(z);
% %         
% %         
% %         
% %         [U,earray1] = xupdateCG(b,param,z,opts,U,Lam, 1e-7,10);
% %         
% %         e = param.E*U - b;
% %         
% %         e_ft = itemp_FT(U);
% %         
% %         
% %         cost = [cost, sum(abs(e(:)).^2)  +  sum(abs(e_ft(:)))*opts.mu];
% %         
% %         if in>1
% %             if abs(cost(end) - cost(end-1))/abs(cost(end)) < 1e-3
% %                 break;
% %             end
% %         end
% %         
% %         
% %         
% %     end
% %     opts.beta=opts.beta*opts.betarate;
% % end
% % Zhat_ktslr=U;
% % Time_ktslr=toc;
% % NMSE_ktslr=RMSE_modi(U,Xtrue);
% % similarity_index=[];
% % for i =1:1:nt
% %     mssim=ssim(abs(Zhat_ktslr(:,:,i)/max(max(Zhat_ktslr(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
% %     similarity_index(i)=mssim;
% % end
% % sim_ktslr=min(similarity_index);
% % save('C:\Users\sbabu\Desktop\Results_MC1\lowres\Xhat_ktslr','Zhat_ktslr');
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LplusS-Otazo perf param %%%%%%%%%%%%%%%%
% % [nx,ny,nt,nc]=size(kdata);
% % tmp = sqrt(sum(abs((b1)).^2,3));
% % b1c = div0(b1,tmp);
% % tmp = max(Xtrue(:));
% % Xtrue = div0(Xtrue,tmp);
% % param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);
% % % % A = @(z)A_fhp3D_p(z);
% % % % At = @(z) At_fhp3D_p(z);
% % % % param.A=A;
% % % % param.At=At;
% % param.d=kdata;
% % 
% % param.T=TempFFT(3);
% % tic;
% % param.lambda_L=0.01;
% % param.lambda_S=0.01;
% % param.nite=50;
% % param.tol=0.0025;
% % M=param.E'*param.d;
% % M=reshape(M,[nx*ny,nt]);
% % Lpre=M;
% % S=zeros(nx*ny,nt);
% % ite=0;
% % while(1)
% %     ite=ite+1;
% %     M0=M;
% %     [Ut,St,Vt]=svd(M-S,0);
% %     St=diag(SoftThresh(diag(St),St(1)*param.lambda_L));
% %     L=Ut*St*Vt';
% %     S=reshape(param.T'*(SoftThresh(param.T*reshape(M-Lpre,[nx,ny,nt]),param.lambda_S)),[nx*ny,nt]);
% %     resk=param.E*reshape(L+S,[nx,ny,nt])-param.d;
% %     M=L+S-reshape(param.E'*resk,[nx*ny,nt]);
% %     Lpre=L;
% %     tmp2=param.T*reshape(S,[nx,ny,nt]);
% %     if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
% % end
% % Xhat_LpS1=L+S;
% % Xhat_LpS=reshape(Xhat_LpS1,[nx,ny,nt]);
% % Time_Otazo= toc;
% % NMSE_Otazo=RMSE_modi(Xhat_LpS,Xtrue);
% % similarity_index=[];
% % for i =1:1:nt
% %     mssim=ssim(abs(Xhat_LpS(:,:,i)/max(max(Xhat_LpS(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
% %     similarity_index(i)=mssim;
% % end
% % sim_Otazo=min(similarity_index)
% % save('C:\Users\sbabu\Desktop\Results_MC1\lowres\Xhat_LS_otazo','Xhat_LpS');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LplusS-lin%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  tic;
% % [nx,ny,nt,nc]=size(kdata);
% % tmp = sqrt(sum(abs((b1)).^2,3));
% % %  tmp = max(kdata(:));
% % % kdata = div0(kdata,tmp);
% % b1c = div0(b1,tmp);
% %  tmp = max(Xtrue(:));
% %  Xtrue = div0(Xtrue,tmp);
% % param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);
% % 
% % %   param.d=param.E*Xtrue;
% %  param.d=kdata;
% % param.T = getT(nx,ny,nt);
% % param.nite=10;
% % param.scaleL = 130/1.2775;
% % param.scaleS = 1/1.2775;
% % param.lambda_L=0.01;
% % param.lambda_S=0.01*param.scaleS;
% % param.Xinf = zeros(nx*ny,nt);
% % %% POGM
% % [L_pogm,S_pogm,xdiff_pogm,cost_pogm,time_pogm,rankL_pogm] = PGM(param,'pogmS',1,'pogmL',1);
% % %% Display: 4 frames
% % L = L_pogm;S = S_pogm;
% % LplusS=L+S;
% % Time_lin=toc;
% % % NMSE_lin=RMSE_modi(LplusS,Xtrue);
% % % similarity_index=[];
% % % for i =1:1:nt
% % %     mssim=ssim(abs(LplusS(:,:,i)/max(max(LplusS(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
% % %     similarity_index(i)=mssim;
% % % end
% % % sim_lin=min(similarity_index);
% % save('C:\Users\sbabu\Desktop\Results_MC1\lowres\Xhat_LS_lin','LplusS');
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AltGDmin %%%%%%%%%%
 tic;
 [nx,ny,nt,nc]=size(kdata);
% % tmp = sqrt(sum(abs((b1)).^2,3));
% % b1c = div0(b1,tmp);
% % %  tmp = max(kdata(:));
% % % kdata = div0(kdata,tmp);
% % param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);
% % %y=param.E*Xtrue;
% %  y=kdata;
% % [nx,ny,nt,nc]=size(kdata);
% % % mask=kdata(:,:,:,1)~=0;
% % % mask1=repmat(mask,[1,1,1,nc]);
% % 
% % % 
% % % [nx,ny,nt,nc]=size(y);
% % n=nx*ny;
% % [zbar_hat,flag,resNE,iter] = cgls_mean(@E_forw_for_mean,@E_back_for_mean, y,0,1e-36,10);
% % ybar_hat=param.E*repmat(zbar_hat,[1,1,nt]);
% % yinter=y-ybar_hat;
% % sum_mk=nnz(y);
% % for k=1:1:nt
% %     mk(k)=nnz(y(:,:,k,:));
% % end
% % m=max(mk);
% % %%%%%%%%%%%%%%%%%%%%% Init Start %%%%%%%%%%%%%%%%%%
% % C_tilda=36;
% % alpha=C_tilda*norm(yinter(:))^2/sum_mk;
% % Y_trunc=yinter;
% % Y_trunc(abs(yinter)>sqrt(alpha))=0;
% % X0_temp=param.E'*Y_trunc;
% % DiagMat=diag(ones(1,nt)./sqrt(mk*m));
% % X0_image=X0_temp;
% % X0=(reshape(X0_image,[nx*ny,nt]))*(DiagMat);
% % %%%%%%%%%%%%%%%%%%%%%% Automated r %%%%%%%%%%%%%%
% % r_big=floor(min([n/10,nt/10,m/10]));
% % [Utemp, Stemp,~]=svds(X0,r_big);
% % SS=diag(Stemp);
% % E=sum(SS.^2);
% % Esum=0;
% % for i=1:1:r_big
% %     Esum=Esum+((SS(i))^2);
% %     if Esum >(E*0.85)
% %         break
% %     end
% % end
% % 
% % r=i+1;
% % r=min(r,r_big);
% % U0=Utemp(:,1:r);
% % %%%%%%%%%%%%%%%%%%
% % % eta=(Stemp(1)-Stemp(r))/(Stemp(1)+Stemp(r));
% % %%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%Init End %%%%%%%%%%%%%%%%%%%%%%%
% % Uhat=U0;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%% altGDmin Start %%%%%%%%%%%%%%%%%%
% % T=70;
% % %   AU=zeros(nx,ny,nt,nc,r);
% %    y_temp=reshape(yinter,[nx*ny,nt,nc]);
% %    
% % for t = 1 : T
% %     Uhatm=reshape(Uhat,[nx,ny,r]);
% %     B = E_forw_for_AU_new(Uhatm,y_temp);
% %     X=reshape(Uhat*B,[nx,ny,nt]);
% %     Z=param.E'*((param.E*X)-yinter);
% %     Z_mat=reshape(Z,[nx*ny,nt]);
% %     Grad_U=Z_mat*B';
% %     if t==1
% %         eta=1/(7*norm(Grad_U));
% %     end
% %     Uhat_t0=Uhat;
% %     Uhat=Uhat-eta*Grad_U;
% %     [Qu,~]  =  qr(Uhat,0);
% %     Uhat  =  Qu(:, 1:r);
% %     Uhat_t1=Uhat;
% %     Subspace_d= ( norm((Uhat_t0 - Uhat_t1*(Uhat_t1'*Uhat_t0)), 'fro')/sqrt(r));
% %     if  (Subspace_d <=.01)
% %         break;
% %     end
% % end
% % 
% % X_GD=X+zbar_hat;
% % Time_GD=toc;
% NMSE_GD=RMSE_modi(X_GD,Xtrue);
% similarity_index=[];
% for i =1:1:nt
%     mssim=ssim(abs(X_GD(:,:,i)/max(max(X_GD(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
%     similarity_index(i)=mssim;
% end
% sim_GD=min(similarity_index);
% save('C:\Users\sbabu\Desktop\Results_MC\0016\Xhat_GD','X_GD');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AltGDmin MRI %%%%%%%%%%
% param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);
% tic;
% 
% % y=param.E*Xtrue;
%  tmp = max(kdata(:));
% kdata = div0(kdata,tmp);
%  y=kdata;
% [nx,ny,nt,nc]=size(kdata);
% mask=kdata(:,:,:,1)~=0;
% mask1=repmat(mask,[1,1,1,nc]);
% 
% [nx,ny,nt,nc]=size(y);
% n=nx*ny;
% [zbar_hat,flag,resNE,iter] = cgls_mean(@E_forw_for_mean,@E_back_for_mean, y,0,1e-36,10);
% ybar_hat=param.E*repmat(zbar_hat,[1,1,nt]);
% yinter=y-ybar_hat;
% sum_mk=nnz(y);
% for k=1:1:nt
%     mk(k)=nnz(y(:,:,k,:));
% end
% m=max(mk);
% C_tilda=36;
% alpha=C_tilda*norm(yinter(:))^2/sum_mk;
% Y_trunc=yinter;
% Y_trunc(abs(yinter)>sqrt(alpha))=0;
% X0_temp=param.E'*Y_trunc;
% DiagMat=diag(ones(1,nt)./sqrt(mk*m));
% X0_image=X0_temp;
% X0=(reshape(X0_image,[nx*ny,nt]))*(DiagMat);
% r_big=floor(min([n/10,nt/10,m/10]));
% [Utemp, Stemp,~]=svds(X0,r_big);
% SS=diag(Stemp);
% % tic;
% E=sum(SS.^2);
% Esum=0;
% for i=1:1:r_big
%     Esum=Esum+((SS(i))^2);
%     if Esum >(E*0.85)
%         break
%     end
% end
% r=i+1;
% r=min(r,r_big);
% U0=Utemp(:,1:r);
% Uhat=U0;
% T=200;
% %   AU=zeros(nx,ny,nt,nc,r);
%   y_temp=reshape(yinter,[nx*ny,nt,nc]);
% for t = 1 : T
%     Uhatm=reshape(Uhat,[nx,ny,r]);
%     B = E_forw_for_AU_new(Uhatm,y_temp);
%     X=reshape(Uhat*B,[nx,ny,nt]);
%     Z=param.E'*((param.E*X)-yinter);
%     Z_mat=reshape(Z,[nx*ny,nt]);
%     Grad_U=Z_mat*B';
%     if t==1
%         eta=1/(7*norm(Grad_U));
%     end
%     Uhat_t0=Uhat;
%     Uhat=Uhat-eta*Grad_U;
%     [Qu,~]  =  qr(Uhat,0);
%     Uhat  =  Qu(:, 1:r);
%     Uhat_t1=Uhat;
% %     error(t)=norm(Xtrue(:)-X(:))
%       Subspace_d= ( norm((Uhat_t0 - Uhat_t1*(Uhat_t1'*Uhat_t0)), 'fro')/sqrt(r));
% %     if  (Subspace_d <=.01)
% %         break;
% %     end
% end
% X_GD=X+zbar_hat;
% yk=y-(param.E*X_GD);
% global k
% % tic;
% for k=1:nt
%     Ehat(:,:,k)=cgls_modi(@E_forw_for_Ak,@E_back_for_Ak,squeeze(yk(:,:,k,:)),0,1e-36,3);
% end
% % Time=toc;
% % tic
% %  [Ehat,flag,resNE,iter] = cgls_modi(@E_forw_for_mec,@E_back_for_mec,yk,0.1,1e-36,3);
% Time_MRI=toc;
%  X_MRI=X+zbar_hat+Ehat;
% NMSE_MRI=RMSE_modi(X_MRI,Xtrue);
% for i =1:1:nt
%     mssim=ssim(abs(X_MRI(:,:,i)/max(max(X_MRI(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
%     similarity_index(i)=mssim;
% end
% sim_MRI=min(similarity_index);
% 
% save('C:\Users\sbabu\Desktop\Results_MC\0016\Xhat_MRI','X_MRI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AltGDmin-MRI2 %%%%%%%%%%
tic;
tmp = sqrt(sum(abs((b1)).^2,3));
b1c = div0(b1,tmp);
param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);
% y=param.E*Xtrue;
 tmp = max(kdata(:));
kdata = div0(kdata,tmp);
y=kdata;

[nx,ny,nt,nc]=size(kdata);
mask=kdata(:,:,:,1)~=0;
mask1=repmat(mask,[1,1,1,nc]);


[nx,ny,nt,nc]=size(y);
n=nx*ny;
[zbar_hat,flag,resNE,iter] = cgls_mean(@E_forw_for_mean,@E_back_for_mean, y,0,1e-36,10);
ybar_hat=param.E*repmat(zbar_hat,[1,1,nt]);
yinter=y-ybar_hat;
sum_mk=nnz(y);
mk=[];
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
Uhat=U0;
T=70;
% AU=zeros(nx,ny,nt,nc,r);
 y_temp=reshape(yinter,[nx*ny,nt,nc]);
for t = 1 : T
    Uhatm=reshape(Uhat,[nx,ny,r]);
     
    B = E_forw_for_AU_new(Uhatm,y_temp);
% B1=Uhat'*reshape(param.E'*(yinter),[nx*ny,nt]);
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


X_GD=X+zbar_hat;


param.d = y;
param.T=getT(nx,ny,nt);
param.lambda_L=0.01;
param.nite=10;
param.tol=0.0025;
M=param.E'*param.d;
Lpre=M;
Ehat=zeros(nx,ny,nt);
L=X_GD;
param.lambda_S=0.001*max(max(abs(M-L)));
ite=0;
while(1)
    ite=ite+1;
    M0=M;
    Ehat=param.T'*(SoftThresh(param.T*reshape(M-Lpre,[nx,ny,nt]),param.lambda_S));
    resk=param.E*(reshape(L+Ehat,[nx,ny,nt]))-param.d;
    M=L+Ehat-reshape(param.E'*resk,[nx,ny,nt]);
    Lpre=L;
    tmp2=param.T*reshape(Ehat,[nx,ny,nt]);
    if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
end

Time_MRI2=toc;
 X_MRI2=X_GD+Ehat;
% NMSE_MRI2=RMSE_modi(X_MRI2,Xtrue);
% for i =1:1:nt
%     mssim=ssim(abs(X_MRI2(:,:,i)/max(max(X_MRI2(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
%     similarity_index(i)=mssim;
% end
% sim_MRI2=min(similarity_index);
save('C:\Users\sbabu\Desktop\Results_MC1\lowres\Xhat_MRI2','X_MRI2');
fprintf(fid, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f) & %8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial,NMSE_ktslr,Time_ktslr,NMSE_Otazo,Time_Otazo,NMSE_lin, Time_lin,NMSE_MRI2,Time_MRI2);
fprintf(fid2, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)  \n', name, radial,sim_ktslr,Time_ktslr,sim_Otazo,Time_Otazo,sim_lin,Time_lin,sim_MRI2,Time_MRI2);
end
% function AUmat = E_forw_for_AU(U_im)
% global nx ny nt nc r b1c samp
% s = zeros(nx,ny,nt,nc,r); 
% 
% U_im_rep = repmat(U_im, [1 1 1 nt]); %nx ny r q 
% 
% s = bsxfun(@times,U_im_rep,reshape(b1c,[nx,ny,1,1,nc]));
% 
% S=fft2c_mri(s); %nx ny r nt nc
% 
% Spermute = permute(S,[1, 2, 4, 5,3]);% nx ny nt r nc
% 
% samp_r = samp ;% nx ny nt 
% 
% if ~isempty(samp) %cl: samp mask only when cartesian samp
%     %S = bsxfun(@times,S,arg.samp);
%     AUmat = bsxfun(@times,Spermute,samp_r);
% end
% end

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
% function AUmat = E_forw_for_AUk(U_im)
% global nx ny  nc r b1c samp k
% s = zeros(nx,ny,nc,r); 
% 
% % U_im_rep = repmat(U_im, [1 1 1 nt]); %nx ny r q 
% 
% s = bsxfun(@times,U_im,reshape(b1c,[nx,ny,1,nc]));
% 
% S=fft2c_mri(s); %nx ny r nt nc
% 
% Spermute = permute(S,[1, 2, 4,3]);% nx ny nc r
% 
% samp_k = squeeze(samp(:,:,k)) ;% nx ny 
% 
% if ~isempty(samp) %cl: samp mask only when cartesian samp
%     %S = bsxfun(@times,S,arg.samp);
%     AUmat = bsxfun(@times,S,reshape(samp_k,[nx,ny,1,1]));
% end
% % 'if above does not work then repmat arg.samp nt x nc times followed by .* with Spermute' 
% end



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


function Amat_z  = E_forw_for_mec(z) %z: nx ny nt
global nx ny nt nc  b1c samp 

smaps = b1c; % nx ny nc
samp = reshape(samp, [nx,ny,nt,1]) ; %nx ny nt 1
s = bsxfun(@times, z, reshape(b1c,[nx,ny,1,nc]));  % nx ny nt nc
S=fft2c_mri(s); %nx ny nt nc
% Spermute = reshape(S, [nx,ny,1,nc]);  
Amat_z = bsxfun(@times,S,samp);  % nx ny nt nc
end

 
function z  = E_back_for_mec(Y_in) %Y_in: nx ny nt nc
global nx ny nt nc  b1c samp 
smaps = b1c; % nx ny nc
samp = reshape(samp, [nx,ny,nt,1]) ; %nx ny nt 1

S = bsxfun(@times,Y_in,samp); %nx ny nt nc
s = ifft2c_mri(S);  %nx ny nt nc
z = sum(bsxfun(@times,s,reshape(conj(smaps),[nx,ny,1,nc])),4);  %nx ny nt after the sum step 
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