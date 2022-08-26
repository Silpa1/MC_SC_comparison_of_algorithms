
clc;clear all;close all;
global nx ny nt nc  b1c samp mk
[fid,msg] = fopen('Comparison_error.txt','wt');
fprintf(fid, '%s(%s) & %s & %s &  %s & %s & %s   \n','Dataset','Radial','ktslr','L+S-Otazo','L+S-Lin','altGDmin-MRI','altGDmin-MRI2');

% Loading fully sampled kspace data nx ny nt nc and coil sensitivity b1c nx
% ny nc
%kdata full-sampled nc-coil kspace data of size nx x ny x nt x nc   ; b1c: coil sensitivity of size nx x ny x nc 
% kdata_u: simulated pseudo-radially undersampled kspace data of size nx x ny x nt x nc  with zeros where no samples was taken.
load('brain_fullsampled_kspace.mat')
[nx,ny,nt,nc]=size(kdata);
%The true image
name='brain';
im = sqrt(nx*ny)*ifft2(permute(kdata,[1,2,4,3]));
Xtrue = squeeze(sum(im.*repmat(conj(b1c),[1 1 1 nt]),3));
radial=[4,8,16];
for ii=1:1:length(radial)
    % Undersampled kspace Data 
    samp = goldencart(nx,ny,nt,radial(ii));
    mask=repmat(samp,[1,1,1,nc]);
     kdata_u=E_forw(Xtrue);
     % Calling the Algorithms
    [Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata_u,b1c,samp);
    NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
    [Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata_u,b1c,samp);
    NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
    [Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata_u,b1c,samp);
    NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
    [Xhat_GD_mean,Time_GD_mean]=mean_altGDmin_func(kdata_u,b1c,samp);
    NMSE_GD_mean=RMSE_modi(Xhat_GD_mean,Xtrue);
    [Xhat_mri,Time_mri]=altGDmin_mri_func(kdata_u,b1c,samp);
    NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
    [Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata_u,b1c,samp);
    NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
    fprintf(fid, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f) & %8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial,NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin, Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading fully sampled kspace data nx ny nt nc and coil sensitivity b1c nx
% ny nc
load('lowres_fullsampled_kspace.mat')
[nx,ny,nt,nc]=size(kdata_u);
%The true image
im = sqrt(nx*ny)*ifft2(permute(kdata_u,[1,2,4,3]));
Xtrue = squeeze(sum(im.*repmat(conj(b1c),[1 1 1 nt]),3));
radial=[4,8,16];
for ii=1:1:length(radial)
    % Undersampled kspace Data 
    samp = goldencart(nx,ny,nt,radial(ii));
    mask=repmat(samp,[1,1,1,nc]);
    param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
     kdata_u=param.E*Xtrue;
     % Calling the Algorithms
    [Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata_u,b1c,samp);
    NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
    [Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata_u,b1c,samp);
    NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
    [Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata_u,b1c,samp);
    NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
    [Xhat_GD_mean,Time_GD_mean]=mean_altGDmin_func(kdata_u,b1c,samp);
    NMSE_GD_mean=RMSE_modi(Xhat_GD_mean,Xtrue);
    [Xhat_mri,Time_mri]=altGDmin_mri_func(kdata_u,b1c,samp);
    NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
    [Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata_u,b1c,samp);
    NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
    fprintf(fid, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f) & %8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial,NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin, Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);
end
% filenames={ 'multi_coil_lowres_speech.mat'}
% load(filenames{1})
% [~,name,~] = fileparts(filenames{1});
% no_comp=8;
% [k]= coil_compress_withpca(k,no_comp);
% [nx,ny,nc,nt]=size(k);
% im = sqrt(nx*ny)*ifft2(k);
% xm = mean(im,4);
% csm = ismrm_estimate_csm_walsh_modified(xm);
% b1=csm;
% tmp = sqrt(sum(abs((b1)).^2,3));
% b1c = div0(b1,tmp);
% Xtrue = squeeze(sum(im.*repmat(conj(csm),[1 1 1 nt]),3));
% radial=[4,8,16];
% for ii=1:1:length(radial)
%     samp = goldencart(nx,ny,nt,radial(ii));
%     mask=repmat(samp,[1,1,1,nc]);
%     param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
%     kdata=param.E*Xtrue;
%     [Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata,b1c,samp);
%     NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
%     [Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata,b1c,samp);
%     NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
%     [Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata,b1c,samp);
%     NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
%     [Xhat_GD_mean,Time_GD_mean]=mean_altGDmin_func(kdata,b1c,samp);
%     NMSE_GD_mean=RMSE_modi(Xhat_GD_mean,Xtrue);
%     [Xhat_mri,Time_mri]=altGDmin_mri_func(kdata,b1c,samp);
%     NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
%     [Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata,b1c,samp);
%     NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
%     fprintf(fid, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f) & %8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial,NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin, Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filenames={ 'multicoil_ungated_cmr_perf.mat'}
% load(filenames{1})
% [~,name,~] = fileparts(filenames{1});
% no_comp=8;
% [kSpace]= coil_compress_withpca(kSpace,no_comp);
% [nx,ny,nc,nt]=size(kSpace);
% im = sqrt(nx*ny)*ifft2(kSpace);
% xm = mean(im,4);
% csm = ismrm_estimate_csm_walsh_modified(xm);
% b1=csm;
% tmp = sqrt(sum(abs((b1)).^2,3));
% b1c = div0(b1,tmp);
% for i =1:1:nc
%     kdata(:,:,:,i)=kSpace(:,:,i,:);
%  end
% Xtrue = squeeze(sum(im.*repmat(conj(csm),[1 1 1 nt]),3));
% tmp = max(Xtrue(:));
% Xtrue = div0(Xtrue,tmp);
% radial=[4];
% for ii=1:1:length(radial)
%     samp = goldencart(nx,ny,nt,radial(ii));
%     mask=repmat(samp,[1,1,1,nc]);
%     param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
%     kdata=param.E*Xtrue;
%     [Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata,b1c,samp);
%     NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
%     [Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata,b1c,samp);
%     NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
%     [Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata,b1c,samp);
%     NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
%     [Xhat_GD_mean,Time_GD_mean]=mean_altGDmin_func(kdata,b1c,samp);
%     NMSE_GD_mean=RMSE_modi(Xhat_GD_mean,Xtrue);
%     [Xhat_mri,Time_mri]=altGDmin_mri_func(kdata,b1c,samp);
%     NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
%     [Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata,b1c,samp);
%     NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
%     fprintf(fid, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f) & %8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial,NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin, Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filenames={ 'fs_0016_3T_ksp_csm.mat'}
% load(filenames{1})
% [~,name,~] = fileparts(filenames{1});
% k = permute(kspace,[3 4 2 1]);
% [nx,ny,nc,nt]=size(k);
% for c = 1:nc,
%     for t = 1:nt,
%         k(:,:,c,t) = fft2(ifftshift(ifft2(fftshift(k(:,:,c,t)))));
%     end
% end
% no_comp=8;
% [k]= coil_compress_withpca(k,no_comp);
% [nx,ny,nc,nt]=size(k);
% im = sqrt(nx*ny)*ifft2(k);
% xm = mean(im,4);
% csm = ismrm_estimate_csm_walsh_modified(xm);
% b1=csm;
% tmp = sqrt(sum(abs((b1)).^2,3));
% b1c = div0(b1,tmp);
% for i =1:1:nc
%     kdata(:,:,:,i)=k(:,:,i,:);
%  end
% Xtrue= squeeze(sum(im.*repmat(conj(csm),[1 1 1 nt]),3));
% radial=[4];
% for ii=1:1:length(radial)
%     samp = goldencart(nx,ny,nt,radial(ii));
%     mask=repmat(samp,[1,1,1,nc]);
%     param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
%     kdata=param.E*Xtrue;
%     [Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata,b1c,samp);
%     NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
%     [Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata,b1c,samp);
%     NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
%     [Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata,b1c,samp);
%     NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
%     [Xhat_GD_mean,Time_GD_mean]=mean_altGDmin_func(kdata,b1c,samp);
%     NMSE_GD_mean=RMSE_modi(Xhat_GD_mean,Xtrue);
%     [Xhat_mri,Time_mri]=altGDmin_mri_func(kdata,b1c,samp);
%     NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
%     [Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata,b1c,samp);
%     NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
%     fprintf(fid, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f) & %8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial,NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin, Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);
% 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filenames={ 'fs_0019_3T_ksp_csm.mat'}
load(filenames{1})
[~,name,~] = fileparts(filenames{1});
k = permute(kspace,[3 4 2 1]);
[nx,ny,nc,nt]=size(k);
for c = 1:nc,
    for t = 1:nt,
        k(:,:,c,t) = fft2(ifftshift(ifft2(fftshift(k(:,:,c,t)))));
    end
end
no_comp=8;
[k]= coil_compress_withpca(k,no_comp);
[nx,ny,nc,nt]=size(k);
im = sqrt(nx*ny)*ifft2(k);
xm = mean(im,4);
csm = ismrm_estimate_csm_walsh_modified(xm);
b1=csm;
tmp = sqrt(sum(abs((b1)).^2,3));
b1c = div0(b1,tmp);
for i =1:1:nc
    kdata_u(:,:,:,i)=k(:,:,i,:);
 end
Xtrue= squeeze(sum(im.*repmat(conj(csm),[1 1 1 nt]),3));
radial=[4,8,16];
for ii=1:1:length(radial)
    samp = goldencart(nx,ny,nt,radial(ii));
    mask=repmat(samp,[1,1,1,nc]);
    param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
    kdata_u=param.E*Xtrue;
    [Xhat_ktslr,Time_ktslr ]=ktslr_func(kdata_u,b1c,samp);
    NMSE_ktslr=RMSE_modi(Xhat_ktslr,Xtrue);
    [Xhat_LpS_otazo,Time_otazo ]=lpluss_otazo_func(kdata_u,b1c,samp);
    NMSE_otazo=RMSE_modi(Xhat_LpS_otazo,Xtrue);
    [Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata_u,b1c,samp);
    NMSE_lin=RMSE_modi(Xhat_LpS_lin,Xtrue);
    [Xhat_GD_mean,Time_GD_mean]=mean_altGDmin_func(kdata_u,b1c,samp);
    NMSE_GD_mean=RMSE_modi(Xhat_GD_mean,Xtrue);
    [Xhat_mri,Time_mri]=altGDmin_mri_func(kdata_u,b1c,samp);
    NMSE_mri=RMSE_modi(Xhat_mri,Xtrue);
    [Xhat_mri2,Time_mri2]=altGDmin_mri2_func(kdata_u,b1c,samp);
    NMSE_mri2=RMSE_modi(Xhat_mri2,Xtrue);
    fprintf(fid, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f) & %8.4f (%5.2f)& %8.4f (%5.2f) \n', name, radial,NMSE_ktslr,Time_ktslr,NMSE_otazo,Time_otazo,NMSE_lin, Time_lin,NMSE_mri,Time_mri,NMSE_mri2,Time_mri2);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



