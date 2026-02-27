function [V_MB,b]=mOPG(X,YY,h,dimx)
[n,m]=size(X);
[yn,ym]=size(YY);
M=zeros(m,m);
D=m;
%===========step1=======================
for kk=1:ym
    Y=YY(:,kk);
    for j=1:n
        pp=0;s11=zeros(D+1,D+1);s12=zeros(D+1,1);
        for im11=1:n
            u1=X-repmat(X(j,:),n,1);
            v1=sqrt(sum((u1.*u1)'))/h;
            Ker10=1/(sqrt(2*pi))*exp(-1/2*v1.^2);
            Ker1_fm=1/(n*h^m)*sum(Ker10);
            
            
            u2=X(im11,:)-X(j,:);
            v2=sqrt(u2*u2')/h;
            Ker_fz=1/(h^m)*1/(sqrt(2*pi))*exp(-1/2*v2^2);
            wij(im11,j)=Ker_fz/Ker1_fm;%后续用的权重wij
            Wij=Ker_fz/Ker1_fm;
            s11=s11+Wij*([1 u2]'*[1 u2]);
            pp=pp+1/(n*h^m)*1/(sqrt(2*pi))*exp(-1/2*v2^2);
            
            
                s12=s12+Wij*[1 u2]'*Y(im11);
            
        end
        s1=s11;
        s2=s12;
        ab=s1^(-1)*s2;
        aj=ab(1,1);
        bj=ab(2:D+1,:);
        
        M=M+(bj*bj');
    end
end
[b,eigval_MB]=eig(M);
[V_MB,P_MB]=sort(diag(abs(eigval_MB)),'descend');
b=b(:,P_MB);
b=b(:,1:dimx);
for jj=1:dimx
    b(:,jj)= b(:,jj)/sqrt(sum(b(:,jj).^2));
end
end