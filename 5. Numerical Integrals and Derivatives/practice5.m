syms x1
f=1/x1^2*sin(2*pi/x1);

%Gauss-Legendre
P(1)=1+0*x1;
for n=1:5
    if n==1
        P(n+1)=(2*n-1)*x1*P(n)/(n);
    else
        P(n+1)=((2*n-1)*x1*P(n)-(n-1)*P(n-1))/(n);
    end
    P(n+1)=collect(P(n+1));
end
for j=1:4
    x(:,j)=(sort(double(solve(P(6)==0))))/4+0.75+j*0.5;
end
x(:,5)=sort(double(solve(P(6)==0)));
 for n=1:5
        A(n)=double(subs(2/5/(P(5)*diff(P(6),x1)),x1,x(n,5)))/4;
 end
 for j=1:4
    for n=1:5
        fk(n,j)=double(subs(f,x1,x(n,j)));
    end
    I1(j)=dot(A,fk(:,j));
end
IGL=sum(I1);

%Romberg
%for j=0:7
%    h(j+1)=2^(-j)*(3-1);
%    for l=1:2^j
%        fh(l,j+1)=double(subs(f,x1,1+(l-1/2)*h(j+1)));
%    end
%end
eps=1e-7;
T(1,1)=double(subs(f,x1,1)+subs(f,x1,3));
fh(1,1)=0;
k=2;
while k>=2
    
    h(k)=2^(-k+1)*(3-1);
    for l=1:2^(k-1)
        fh(l,k)=double(subs(f,x1,1+(l-1/2)*h(k)));
    end
    
    T(k,1)=1/2*(T(k-1,1)+h(k-1)*sum(fh(1:2^(k-1),k-1)));
    for j=1:k-1
        T(k-j,1+j)=(4^j*T(k-j+1,j)-T(k-j,j))/(4^j-1);
    end
    if abs(T(1,k)-T(1,k-1))>=eps
        k=k+1;
    else
        break
    end
    
end
IR=T(1,k);

Ref=double(int(f,x1,1,3));
errGL=IGL-Ref;
errR=IR-Ref;