function  [Xhat_LpS_lin,Time_lin]=lpluss_lin_func(kdata,b1c,samp)
 tic;
[nx,ny,nt,nc]=size(kdata);
param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);
param.d=kdata;
param.T = getT(nx,ny,nt);
param.nite=10;
param.scaleL = 130/1.2775;
param.scaleS = 1/1.2775;
param.lambda_L=0.01;
param.lambda_S=0.01*param.scaleS;
param.Xinf = zeros(nx*ny,nt);
[L_pogm,S_pogm,xdiff_pogm,cost_pogm,time_pogm,rankL_pogm] = PGM(param,'pogmS',1,'pogmL',1);
L = L_pogm;S = S_pogm;
Xhat_LpS_lin=L+S;
Time_lin=toc;
end