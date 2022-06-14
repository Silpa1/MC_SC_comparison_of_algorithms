clc;
clear all;close all;
global n1 n2 n mk m  S2 q  kk;

filenames={'brain_T2_T1.mat'};


[fid,msg] = fopen('Comparison_error.txt','wt');
[fid2,msg] = fopen('Comparison_sim.txt','wt');
fprintf(fid, '%s(%s) & %s & %s &  %s &  %s &  %s  &  %s  \n','Dataset','Radial','ktslr','L+S-Otazo','L+S-Lin','altGDmin','altGDmin-MRI','altGDmin-MRI2');
fprintf(fid2, '%s(%s) & %s & %s &  %s &  %s &  %s &  %s   \n','Dataset','Radial','ktslr','L+S-Otazo','L+S-Lin','altGDmin','altGDmin-MRI','altGDmin-MRI2');

for jj = 1:1:numel(filenames)
    S = load(filenames{jj});
    [~,name,~] = fileparts(filenames{jj});
    radial=[4,8,16];
    X_image=double(cell2mat(struct2cell(S)));
    
    save('C:\Users\sbabu\Desktop\Results\low_res_8\X_input.mat', 'X_image');
    [n1,n2,q]=size(X_image);
    n=n1*n2;
    X_mat=reshape(X_image,[n,q]);
    
    for ii=1:1:length(radial)
        GD_MLS_time=0;
        GD_MLS_error=0;
        [mask1]=goldencart(n1,n2,q,radial(ii));
        mask = fftshift(fftshift(mask1,1),2);
        Samp_loc=double(find(logical(mask)));
        mask3=reshape(mask,[n, q]);
        mk=[];
        for i=1:1:q
            mk(i)=length(find(logical(mask3(:,i))));
            S2(1:mk(i),i)=double(find(logical(mask3(:,i))));
        end
        m=max(mk);
        Y=zeros(m,q);
        for k=1:1:q
            ksc = reshape( fft2( reshape(X_mat(:,k), [n1 n2]) ), [n,1]) ;
            Y(1:mk(k),k)=double(ksc(S2(1:mk(k),k)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ktslr %%%%%%%%%%%%%%%%%%%%%%%%%%
        S=find(mask~=0);
        Akt = @(z)A_fhp3D(z,S,n1,n2,q);
        Aktt=@(z)At_fhp3D(z,S,n1,n2,q);
        step_size = [1,1,1];
        [D,Dt] = defDDt(step_size);
        b = Akt(X_image);
        tic;
        x_init = Aktt(b);
        mu1 =1e-10;
        mu2 =4e-9;
        opts.mu1 = mu1;
        opts.mu2 = mu2;
        opts.p=0.1;
        [~,sq,~]=givefastSVD(reshape(x_init, n1*n2,q));
        opts.beta1=10./max(sq(:));
        opts.beta2=10./max(abs(x_init(:)));
        opts.beta1rate = 50;
        opts.beta2rate = 25;
        opts.outer_iter =15;
        opts.inner_iter = 50;
        [Xhat_ktslr,cost,opts] = minSNandTV(Akt,Aktt,D,Dt,x_init,b,1,opts);
        Time_Ktslr=toc;
        Error_Ktslr= RMSE_modi(Xhat_ktslr,X_image);
        similarity_index=[];
        for i =1:1:q
            similarity_index(i)=ssim(abs(Xhat_ktslr(:,:,i)/max(max(Xhat_ktslr(:,:,i)))),abs(X_image(:,:,i)/max(max(X_image(:,:,i)))));
        end
        sim_Ktslr=min(similarity_index)
        %save('C:\Users\sbabu\Desktop\Result\brain_8\Xhat_ktslr.mat', 'Xhat_ktslr');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L+S-Otazo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        param.Samp_loc=Samp_loc;
        A = @(z)A_fhp3D(z, Samp_loc,n1,n2,q);
        At = @(z) At_fhp3D(z, Samp_loc, n1,n2,q);
        param.A=A;
        param.At=At;
        param.d = A(X_mat);
        param.T=TempFFT(3);
        tic;
        param.lambda_L=0.01;
        param.lambda_S=0.01;
        param.nite=50;
        param.tol=0.0025;
        M=At(param.d);
        
        M=reshape(M,[n1*n2,q]);
        Lpre=M;
        S=zeros(n1*n2,q);
        ite=0;
        while(1)
            ite=ite+1;
            M0=M;
            [Ut,St,Vt]=svd(M-S,0);
            St=diag(SoftThresh(diag(St),St(1)*param.lambda_L));
            L=Ut*St*Vt';
            S=reshape(param.T'*(SoftThresh(param.T*reshape(M-Lpre,[n1,n2,q]),param.lambda_S)),[n1*n2,q]);
            resk=param.A(reshape(L+S,[n1,n2,q]))-param.d;
            M=L+S-reshape(param.At(resk),[n1*n2,q]);
            Lpre=L;
            tmp2=param.T*reshape(S,[n1,n2,q]);
            if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
        end
        Xhat_LpS1=L+S;
        Xhat_LpS=reshape(Xhat_LpS1,[n1,n2,q]);
        Time_LSparse= toc;
        Error_LSparse=RMSE_modi(Xhat_LpS,X_image);
        similarity_index=[];
        for i =1:1:q
            mssim=ssim(abs(Xhat_LpS(:,:,i)/max(max(Xhat_LpS(:,:,i)))),abs(X_image(:,:,i)/max(max(X_image(:,:,i)))));
            similarity_index(i)=mssim;
        end
        sim_LpS=min(similarity_index)
        %save('C:\Users\sbabu\Desktop\Result\brain_8\Xhat_LpS.mat', 'Xhat_LpS');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L+S-Lin %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp = max(X_image(:));
        Xtrue = div0(X_image,tmp);
        smap = ones(n1,n2,1);
        
        E=getE_sc(smap,q,'samp',mask1);
        tic;
        d = E*Xtrue;
        % add noise
        rng(0)
        %        dn = randn(size(d)) + 1i * randn(size(d));
        %         param.snr_db =23;
        %         param.scale_noise = norm(d(:)) / norm(dn(:)) / 10.^(param.snr_db / 20);
        param.d = d;
        opt.d = param.d;
        % prepare for regularization scaling
        L = E'*param.d;
        res = E*L-param.d;
        [~,St,~]=svd(reshape(L,[n,q])-reshape(E'*res,[n,q]),0);
        
        param.E = getE(smap,q,'samp',mask1);
        param.T = getT(n1,n2,q);
        param.nite = 10;
        param.scaleL = St(1);
        param.scaleS = 1/1.887;
        param.lambda_L=0.01;
        param.lambda_S=0.05*param.scaleS;
        %param.Xinf = reshape(Xinf.pincat,nx*ny,nt);
        %% ISTA
        %[L_ista,S_ista,xdiff_ista,cost_ista,time_ista,rankL_ista] = PGM(param);
        %% FISTA
        %[L_fista,S_fista,xdiff_fista,cost_fista,time_fista,rankL_fista] = PGM(param,'fistaL',1,'fistaS',1);
        %% POGM
        [L_pogm,S_pogm,xdiff_pogm,cost_pogm,time_pogm,rankL_pogm] = PGM(param,'pogmS',1,'pogmL',1);
        %% Display: 1 frame
        L = L_pogm;S = S_pogm;
        LplusS=L+S;
        Time_LplusS_jeff=toc;
        Error_LplusS_jeff=RMSE_modi(LplusS,Xtrue);
        similarity_index=[];
        for i =1:1:q
            mssim=ssim(abs(LplusS(:,:,i)/max(max(LplusS(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
            similarity_index(i)=mssim;
        end
        %save('C:\Users\sbabu\Desktop\Result\brain_8\Xhat_LpS_lin.mat', 'LplusS');
        sim_LplusS_jeff=min(similarity_index);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%altGDmin%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m=max(mk);
        Y=Y(1:m,1:q);
        T=70;
        tic;
        [Uhat]=initAltGDMin(Y);
        [Uhat2, Bhat2]=GDMin_wi(T,Uhat,Y);
        Xhat_GD=reshape(Uhat2*Bhat2,[n1, n2,q]);
        Time_GD=  toc;
        Error_GD=RMSE_modi(Xhat_GD,X_image);
        similarity_index=[];
        for i =1:1:q
            mssim=ssim(abs(Xhat_GD(:,:,i)/max(max(Xhat_GD(:,:,i)))),abs(X_image(:,:,i)/max(max(X_image(:,:,i)))));
            similarity_index(i)=mssim;
        end
        save('C:\Users\sbabu\Desktop\Result\brain_8\Xhat_GD.mat', 'Xhat_GD');
        sim_GD=min(similarity_index);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AltgdMin + Sparse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T=70;
        tic;
        L=[];
        [Xbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,10);
        Ybar_hat=Afft(Xbar_hat);
        Ybar_hat=reshape(Ybar_hat,[m,q]);
        Yinter=Y-Ybar_hat;
        [Uhat]=initAltGDMin(Yinter);
        [Uhat2, Bhat2]=GDMin_wi(T,Uhat,Yinter);
        xT=Uhat2*Bhat2;
        L(:,1:q)=xT+Xbar_hat;
        
        param.Samp_loc=Samp_loc;
        A = @(z)A_fhp3D(z, Samp_loc,n1,n2,q);
        At = @(z) At_fhp3D(z, Samp_loc, n1,n2,q);
        param.A=A;
        param.At=At;
        param.d = A(X_image(:,:,1:q));
        param.T=TempFFT(3);
        param.lambda_L=0.01;
        
        param.nite=10;
        param.tol=0.0025;
        M=At(param.d);
        M=reshape(M,[n1*n2,q]);
        Lpre=M;
        S=zeros(n1*n2,q);
        param.lambda_S=0.001*max(max(abs(M-L)));
        ite=0;
        while(1)
            ite=ite+1;
            M0=M;
            % sparse update
            S=reshape(param.T'*(SoftThresh(param.T*reshape(M-Lpre,[n1,n2,q]),param.lambda_S)),[n1*n2,q]);
            % data consistency
            resk=param.A(reshape(L+S,[n1,n2,q]))-param.d;
            M=L+S-reshape(param.At(resk),[n1*n2,q]);
            Lpre=L;
            tmp2=param.T*reshape(S,[n1,n2,q]);
            if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
        end
        Xhat_MGDS1=L+S;
        Xhat_MGDS=reshape(Xhat_MGDS1,n1,n2,q);
        Time_GD_Sparse=  toc;
        %save('C:\Users\sbabu\Desktop\Result\brain_8\Xhat_MGDS.mat', 'Xhat_MGDS');
        % Time_GD_Sparse= [Time_GD_Sparse, toc];
        Error_GD_Sparse=RMSE_modi(Xhat_MGDS,X_image);
        similarity_index=[];
        for i =1:1:q
            mssim=ssim(abs(Xhat_MGDS(:,:,i)/max(max(Xhat_MGDS(:,:,i)))),abs(X_image(:,:,i)/max(max(X_image(:,:,i)))));
            similarity_index(i)=mssim;
        end
        sim_MGDS=min(similarity_index)
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AltgdMin + MEC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        m=max(mk);
        L=[];
        T=70;
        tic;
        [Xbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,10);
        Ybar_hat=Afft(Xbar_hat);
        Ybar_hat=reshape(Ybar_hat,[m,q]);
        Yinter=Y-Ybar_hat;
        [Uhat]=initAltGDMin(Yinter);
        [Uhat2, Bhat2]=GDMin_wi(T,Uhat,Yinter);
        xT=Uhat2*Bhat2;
        L(:,1:q)=xT+Xbar_hat;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ymec=Y-Afft(L);
        E_mec=[];
        for kk=1:1:q
            E_mec(:,kk)=cgls_modi(@Afft_modi,@At_modi, Ymec(:,kk) ,0,1e-36,3);
        end
        Xhat_GD_MEC1=L+E_mec;
        Xhat_GD_MEC=reshape(Xhat_GD_MEC1,[n1, n2,q]);
        
        Time_GD_MEC=  toc;
        Error_GD_MEC=RMSE_modi(Xhat_GD_MEC,X_image);
        similarity_index=[];
        for i =1:1:q
            mssim=ssim(abs(Xhat_GD_MEC(:,:,i)/max(max(Xhat_GD_MEC(:,:,i)))),abs(X_image(:,:,i)/max(max(X_image(:,:,i)))));
            similarity_index(i)=mssim;
        end
        sim_GD_MEC=min(similarity_index)
        %save('C:\Users\sbabu\Desktop\Result\brain_8\Xhat_MGD_MEC.mat', 'Xhat_GD_MEC');
        fprintf(fid, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f) & %8.4f (%5.2f) \n', name, radial(ii),Error_Ktslr,Time_Ktslr,Error_LSparse,Time_LSparse,Error_LplusS_jeff,Time_LplusS_jeff,Error_GD,Time_GD,Error_GD_MEC,Time_GD_MEC,Error_GD_Sparse,Time_GD_Sparse);
        fprintf(fid2, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)  \n', name, radial(ii),sim_Ktslr,Time_Ktslr,sim_LpS,Time_LSparse,sim_LplusS_jeff,Time_LplusS_jeff,sim_GD,Time_GD,sim_GD_MEC,Time_GD_MEC,sim_MGDS,Time_GD_Sparse);
 
    end
end
fclose(fid);
fclose(fid2);

function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end
