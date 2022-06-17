clc;
clear all;close all;
global n1 n2 n mk m  S2 q  kk;
filenames=  {'Pincat.mat','brain_T2_T1.mat','speech_seq.mat','Cardiac_ocmr_data.mat','lowres_speech.mat','FB_ungated.mat'};
%filenames={'brain_T2_T1.mat'};


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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L+S-Lin %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp = max(X_image(:));
        Xtrue = div0(X_image,tmp);
        smap = ones(n1,n2,1);
        E=getE_sc(smap,q,'samp',mask1);
        tic;
        d = E*Xtrue;
        param.d = d;
        opt.d = param.d;
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
        
        [L_pogm,S_pogm,xdiff_pogm,cost_pogm,time_pogm,rankL_pogm] = PGM(param,'pogmS',1,'pogmL',1);
        L = L_pogm;S = S_pogm;
        LplusS=L+S;
        Time_LplusS_jeff=toc;
        Error_LplusS_jeff=RMSE_modi(LplusS,Xtrue);
        similarity_index=[];
        for i =1:1:q
            mssim=ssim(abs(LplusS(:,:,i)/max(max(LplusS(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
            similarity_index(i)=mssim;
        end
        sim_LplusS_jeff=min(similarity_index);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%altGDmin%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m=max(mk);
        Y=Y(1:m,1:q);
        T=70;
        tic;
        [Uhat]=initAltGDMin(Y);
        [Uhat2, Bhat2]=AltGDmin(T,Uhat,Y);
        X_hat=reshape(Uhat2*Bhat2,[n1, n2,q]);
        Time_GD=  toc;
        Error_GD=RMSE_modi(X_hat,X_image);
        similarity_index=[];
        for i =1:1:q
            mssim=ssim(abs(X_hat(:,:,i)/max(max(X_hat(:,:,i)))),abs(X_image(:,:,i)/max(max(X_image(:,:,i)))));
            similarity_index(i)=mssim;
        end
        sim_GD=min(similarity_index);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AltgdMin + Sparse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T=70;
        tic;
        L=[];
        [zbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,10);
        Ytemp=reshape(Afft(zbar_hat),[m,q]);
        Ybar=Y-Ytemp;
        [Uhat]=initAltGDMin(Ybar);
        [Uhat2, Bhat2]=AltGDmin(T,Uhat,Ybar);
        X_hat=Uhat2*Bhat2;
       
        
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
        Ehat=zeros(n1*n2,q);
        L(:,1:q)=X_hat+zbar_hat;
        param.lambda_S=0.001*max(max(abs(M-L)));
        ite=0;
       
        while(1)
            ite=ite+1;
            M0=M;
            Ehat=reshape(param.T'*(SoftThresh(param.T*reshape(M-Lpre,[n1,n2,q]),param.lambda_S)),[n1*n2,q]);
            resk=param.A(reshape(L+Ehat,[n1,n2,q]))-param.d;
            M=L+Ehat-reshape(param.At(resk),[n1*n2,q]);
            Lpre=L;
            tmp2=param.T*reshape(Ehat,[n1,n2,q]);
            if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
        end
        Zhat=L+Ehat;
        Zhat_MRI2=reshape(Zhat,n1,n2,q);
        Time_MRI2=  toc;
        Error_MRI2=RMSE_modi(Zhat_MRI2,X_image);
        similarity_index=[];
        for i =1:1:q
            mssim=ssim(abs(Zhat_MRI2(:,:,i)/max(max(Zhat_MRI2(:,:,i)))),abs(X_image(:,:,i)/max(max(X_image(:,:,i)))));
            similarity_index(i)=mssim;
        end
        sim_MRI2=min(similarity_index)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AltgdMin + MEC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m=max(mk);
        L=[];
        T=70;
        tic;
        
        
        [zbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,10);
        Ytemp=reshape(Afft(zbar_hat),[m,q]);
        Ybar=Y-Ytemp;
        [Uhat]=initAltGDMin(Ybar);
        [Uhat2, Bhat2]=AltGDmin(T,Uhat,Ybar);
        X_hat=Uhat2*Bhat2;
        
        Yhat_hat=Y-Afft(X_hat+zbar_hat);
        Ehat=[];
        for kk=1:1:q
            Ehat(:,kk)=cgls_modi(@Afft_modi,@At_modi, Yhat_hat(:,kk) ,0,1e-36,3);
        end
        Zhat=X_hat+zbar_hat+Ehat;
        Zhat_MRI=reshape(Zhat,[n1, n2,q]);
        
        Time_MRI=  toc;
        Error_MRI=RMSE_modi(Zhat_MRI,X_image);
        similarity_index=[];
        for i =1:1:q
            mssim=ssim(abs(Zhat_MRI(:,:,i)/max(max(Zhat_MRI(:,:,i)))),abs(X_image(:,:,i)/max(max(X_image(:,:,i)))));
            similarity_index(i)=mssim;
        end
        sim_MRI=min(similarity_index)
        fprintf(fid, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f) & %8.4f (%5.2f) \n', name, radial(ii),Error_Ktslr,Time_Ktslr,Error_LSparse,Time_LSparse,Error_LplusS_jeff,Time_LplusS_jeff,Error_GD,Time_GD,Error_MRI,Time_MRI,Error_MRI2,Time_MRI2);
        fprintf(fid2, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)  \n', name, radial(ii),sim_Ktslr,Time_Ktslr,sim_LpS,Time_LSparse,sim_LplusS_jeff,Time_LplusS_jeff,sim_GD,Time_GD,sim_MRI,Time_MRI,sim_MRI2,Time_MRI2);
        
    end
end
fclose(fid);
fclose(fid2);

function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end
