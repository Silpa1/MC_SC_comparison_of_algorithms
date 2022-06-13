clc;
clear all;close all;
global n1 n2 n mk m  S2 gg q  jj n3 kk pp ww;
%filenames=  {'brain_T2_T1.mat'} ;
%filenames=  {'CSM.mat'} ;
% filenames=  {'Cardiac_ocmr_data.mat'}%','brain_T2_T1.mat','speech_seq.mat','Cardiac_ocmr_data.mat','FB_ungated.mat'}%,'FB_ungated.mat',,'image_xyt'}  {'Cardiac_ocmr_data.mat'}%; %Cardiac'FB_ungated.mat',

%filenames={'FB_ungated.mat'}%'lowres_speech.mat'}%
%filenames=  {'freebreathing_ungated_Cardiac_cine_Cartesian.mat'} ;
%filenames={'Pincat.mat'};
%filenames={'speech_seq.mat'}%,'lowres_speech.mat'};
%filenames={'image_xyt'};
filenames={'lowres_speech.mat'};
ns=[8];

[fid,msg] = fopen('Subspace_Sparse.txt','wt');
fprintf(fid, '%s & %s    \n','Subspace','AltGDMin MRI');
for jj = 1:1:numel(filenames)
    S = load(filenames{jj});
    X=cell2mat(struct2cell(S));
    [~,name,~] = fileparts(filenames{jj});% Best to load into an output variable.
    radial=[16];
    x=X;
    x=double(x);
    [n1,n2,q]=size(x);
    n3=q;
    X1=reshape(x,[n1*n2,q]);
    n=n1*n2;
    % Table = cell(length(radial)*numel(filenames), 11);
    Error_Ktslr=[];
    Time_Ktslr=[];
    Time_ST=[];
    Error_LSparse=[];
    Error_ST=[];
    Time_GD_Sparse=[];
    
    GD_MLS_time=0;
    GD_MLS_error=0;
    %  [mask] = goldenangle(n1,n2,q,radial(ii));%load('ocmr_test_kspace_mask_modlRecon_16lines.mat','mask_ocmr');
    %[mask] = load('ocmr_test_kspace_mask_modlRecon_16lines.mat','mask_ocmr');
    %mask=cell2mat(struct2cell(mask));%goldencart(n1,n2,q,radial(ii));
    [mask]=goldencart(n1,n2,q,16);
    %[mask]=strucrand(n1,n2,q,radial(ii));
    % mask2=cell2mat(struct2cell(mask2));
    mask = fftshift(fftshift(mask,1),2);
    Samp_loc=double(find(logical(mask)));
    mask3=reshape(mask,[n1*n2, q]);
    mk=[];
    for i=1:1:q
        mk(i)=length(find(logical(mask3(:,i))));
        S2(1:mk(i),i)=double(find(logical(mask3(:,i))));
    end
    m=max(mk);
    Y=zeros(m,q);
    for k=1:1:q
        ksc = reshape( fft2( reshape(X1(:,k), [n1 n2]) ), [n,1]) ;
        Y(1:mk(k),k)=double(ksc(S2(1:mk(k),k)));
    end
    for ii=1:1:length(ns)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subspace Tracking+ AltgdMin + MEC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Mean +LowRank
        %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ subspace^^^^^^^^^^^^^^^^^^^^^^^^^^
        ss=ns(ii);
        n3=q/ss;
        pp=n3;
        Xhat=[];
        Xhat1=[];
        tic;
        for ww=1:1:ss
            gg=mk((ww-1)*n3+1:ww*n3);
            m=max(gg);
            Ys=Y(1:m,(ww-1)*n3+1:ww*n3);
            
            L=[];
            T=5;
            if ww==1
                T=70;
            end
            
            
            
           [Xbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Ys,0,1e-36,10);
           Ybar_hat=Afft(Xbar_hat);
           Ybar_hat=reshape(Ybar_hat,[m,n3]);
            Yinter=Ys-Ybar_hat;
            if ww==1
                [Uhat]=initAltGDMin(Yinter);
                %[Uhat]=initAltGDMin(Ys);
            end
            [Uhat, Bhat]=GDMin_wi(T,Uhat,Yinter);
            %[Uhat, Bhat]=GDMin_wi(T,Uhat,Ys);
            xT=Uhat*Bhat;
            L(:,1:n3)=xT+Xbar_hat;
           % L(:,1:n3)=xT;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Samp_loc=double(find(logical(mask(:,:,(ww-1)*n3+1:ww*n3))));
            param.Samp_loc=Samp_loc;
            A = @(z)A_fhp3D(z, Samp_loc,n1,n2,n3);
            At = @(z) At_fhp3D(z, Samp_loc, n1,n2,n3);
            param.A=A;
            param.At=At;
            param.d = A(X(:,:,(ww-1)*n3+1:ww*n3));
            param.T=TempFFT(3);
            param.lambda_L=0.01;
            param.lambda_S=0.001;
            param.nite=10;
            param.tol=0.0025;
            M=At(param.d);
            M=reshape(M,[n1*n2,n3]);
            Lpre=M;
            S=zeros(n1*n2,n3);
            ite=0;
            while(1)
                ite=ite+1;
                M0=M;
                % sparse update
                S=reshape(param.T'*(SoftThresh(param.T*reshape(M-Lpre,[n1,n2,n3]),param.lambda_S)),[n1*n2,n3]);
                % data consistency
                resk=param.A(reshape(L+S,[n1,n2,n3]))-param.d;
                M=L+S-reshape(param.At(resk),[n1*n2,n3]);
                Lpre=L;
                tmp2=param.T*reshape(S,[n1,n2,n3]);
                if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Xhat=[Xhat,L+S];
            %Xhat1=[Xhat1,S];
        end
        Xhat_mat=reshape(Xhat,[n1, n2,q]);
       % Xhat1_mat=reshape(Xhat1,[n1, n2,q]);
        
        Time_ST=  toc;
        Error_ST=RMSE_modi(Xhat_mat,x);
        %Error_ST=RMSE_modi(Xhat1_mat,x);
        
        %save('C:\Users\sbabu\Desktop\Mini-Batch\Mini_Batch_MEC_Sparse\golden_angle_and_radial_comparison\Xhat_GD_MEC.mat', 'Xhat_GD_MEC');
        similarity_index=[];
        for i =1:1:q
            mssim=ssim(abs(Xhat_mat(:,:,i)/max(max(Xhat_mat(:,:,i)))),abs(x(:,:,i)/max(max(x(:,:,i)))));
            similarity_index(i)=mssim;
        end
        sim_ST=min(similarity_index);
        fprintf(fid, '%d & %8.4f (%5.2f,%8.4f) \n', ss,Error_ST,Time_ST,sim_ST);
    end
end
fclose(fid);




function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end