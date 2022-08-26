function S = E_forw(x)
global nx ny nt nc samp b1c
s = zeros(nx,ny,nt,nc); 
s = bsxfun(@times,x,reshape(b1c,[nx,ny,1,nc]));
S=fft2c_mri(s);
if ~isempty(samp) %cl: samp mask only when cartesian samp
    S = bsxfun(@times,S,samp);
end
end