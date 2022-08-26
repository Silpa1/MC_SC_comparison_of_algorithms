% x = A' * y
function x = E_back(S)
global nx ny nt nc samp b1c
x = zeros(nx,ny,nt);
if ~isempty(samp)
    S = bsxfun(@times,S,samp);
end
s = ifft2c_mri(S);
x = sum(bsxfun(@times,s,reshape(conj(b1c),[nx,ny,1,nc])),4);
end