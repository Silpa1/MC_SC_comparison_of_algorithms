
clc;clear all;close all;
global nx ny nt nc  b1c samp mk
[fid,msg] = fopen('Comparison_error.txt','wt');
[fid2,msg] = fopen('Comparison_sim.txt','wt');
fprintf(fid, '%s(%s) & %s & %s& %s & %s & %s & %s   \n','Dataset','Radial','L+S-ktslr','L+S-Otazo','L+S-Lin','altGDmin_mean','altGDmin-MRI','altGDmin-MRI2');
fprintf(fid2, '%s(%s) & %s & %s& %s & %s & %s & %s   \n','Dataset','Radial','L+S-Lin','altGDmin_mean','altGDmin-MRI','altGDmin-MRI2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kdata=[];
filenames={ 'brain_T2_T1.mat'}
S=load(filenames{1})
[~,name,~] = fileparts(filenames{1});
Xtrue=double(cell2mat(struct2cell(S)));
[nx,ny,nt]=size(Xtrue);
tmp = max(Xtrue(:));
Xtrue = div0(Xtrue,tmp);

radial=[4,8,16];
b1c=ones(nx,ny,1);
for ii=1:1:length(radial)
    samp = goldencart(nx,ny,nt,radial(ii));
    mask=repmat(samp,[1,1,1,nc]);
    param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
    kdata=param.E*Xtrue;
    Single_coil_Algorithms(fid, fid2,name,radial(ii),kdata,b1c,samp,Xtrue)
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% kdata=[];
% filenames={'Cardiac_ocmr_data.mat'}
% S=load(filenames{1})
% [~,name,~] = fileparts(filenames{1});
% Xtrue=double(cell2mat(struct2cell(S)));
% [nx,ny,nt]=size(Xtrue);
% tmp = max(Xtrue(:));
% Xtrue = div0(Xtrue,tmp);
% 
% radial=[4,8,16];
% b1c=ones(nx,ny,1);
% for ii=1:1:length(radial)
%     samp = goldencart(nx,ny,nt,radial(ii));
%     mask=repmat(samp,[1,1,1,nc]);
%     param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
%     kdata=param.E*Xtrue;
%     Single_coil_Algorithms(fid, fid2,name,radial(ii),kdata,b1c,samp,Xtrue)
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kdata=[];
% filenames={ 'speech_seq.mat'}
% S=load(filenames{1})
% [~,name,~] = fileparts(filenames{1});
% Xtrue=double(cell2mat(struct2cell(S)));
% [nx,ny,nt]=size(Xtrue);
% tmp = max(Xtrue(:));
% Xtrue = div0(Xtrue,tmp);
% 
% radial=[4,8,16];
% b1c=ones(nx,ny,1);
% for ii=1:1:length(radial)
%     samp = goldencart(nx,ny,nt,radial(ii));
%     mask=repmat(samp,[1,1,1,nc]);
%     param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
%     kdata=param.E*Xtrue;
%     Single_coil_Algorithms(fid, fid2,name,radial(ii),kdata,b1c,samp,Xtrue)
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kdata=[];
% filenames={ 'FB_ungated.mat'}
% S=load(filenames{1})
% [~,name,~] = fileparts(filenames{1});
% Xtrue=double(cell2mat(struct2cell(S)));
% [nx,ny,nt]=size(Xtrue);
% tmp = max(Xtrue(:));
% Xtrue = div0(Xtrue,tmp);
% 
% radial=[4,8,16];
% b1c=ones(nx,ny,1);
% for ii=1:1:length(radial)
%     samp = goldencart(nx,ny,nt,radial(ii));
%     mask=repmat(samp,[1,1,1,nc]);
%     param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
%     kdata=param.E*Xtrue;
%     Single_coil_Algorithms(fid, fid2,name,radial(ii),kdata,b1c,samp,Xtrue)
% end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kdata=[];
% filenames={ 'lowres_speech.mat'}
% S=load(filenames{1})
% [~,name,~] = fileparts(filenames{1});
% Xtrue=double(cell2mat(struct2cell(S)));
% [nx,ny,nt]=size(Xtrue);
% tmp = max(Xtrue(:));
% Xtrue = div0(Xtrue,tmp);
% 
% radial=[4,8,16];
% b1c=ones(nx,ny,1);
% for ii=1:1:length(radial)
%     samp = goldencart(nx,ny,nt,radial(ii));
%     mask=repmat(samp,[1,1,1,nc]);
%     param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
%     kdata=param.E*Xtrue;
%     Single_coil_Algorithms(fid, fid2,name,radial(ii),kdata,b1c,samp,Xtrue)
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kdata=[];
filenames={ 'Pincat.mat'}
S=load(filenames{1})
[~,name,~] = fileparts(filenames{1});
Xtrue=double(cell2mat(struct2cell(S)));
[nx,ny,nt]=size(Xtrue);
tmp = max(Xtrue(:));
Xtrue = div0(Xtrue,tmp);

radial=[4,8,16];
b1c=ones(nx,ny,1);
for ii=1:1:length(radial)
    samp = goldencart(nx,ny,nt,radial(ii));
    mask=repmat(samp,[1,1,1,nc]);
    param.E=getE(b1c,nt,'samp',mask(:,:,:,1)~=0);
    kdata=param.E*Xtrue;
    Single_coil_Algorithms(fid, fid2,name,radial(ii),kdata,b1c,samp,Xtrue)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
