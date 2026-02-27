function [directions, pgdr_variable] = pgdr(y1,y2,x)

n = size(x,1);
y = [y1 y2];
sigmax=cov(x);
mx=mean(x);
X=(x-repmat(mx,n,1))*sigmax^(-1/2);
del=max(abs(y));
Y=y./(repmat(del,52,1));
h=2;

[eigval_MBx,BX]=mOPG(X,Y,h,size(y,2));
bm = 0.9;s = 0;
K=1;
while s<bm
    s = sum(eigval_MBx(1:K));
    K = K+1;
end
BX=BX(:,1:K);
Dx=K;

for jj=1:10
    [eigval_MBy,BY]=mOPGy(Y,X*BX,h,size(y,2));
    by0=abs(BY);
    BX=swMAVE(X,Y,h,BX,Dx,eigval_MBy,by0);
end

directions = BX;
pgdr_variable = x*BX;

end
