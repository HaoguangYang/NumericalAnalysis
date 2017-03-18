order=3;
step=2000;
syms x                  
Function=x^2*log(2+x);

%Legendre
P(1)=1+0*x;
for n=1:order
    if n==1
        P(n+1)=(2*n-1)*x*P(n)/(n);
    else
        P(n+1)=((2*n-1)*x*P(n)-(n-1)*P(n-1))/(n);
    end
end
for n=1:order+1
    a(n)=int(Function*P(n),x,-1,1)*(2*n-1)/2;
end
s_L=dot(a(:),P(:).');

%Chebyshev
T(1)=1+0*x;
for n=1:order+1
    if n==1
        T(2)=x;
    else
        T(n+1)=2*x*T(n)-T(n-1);%cos(n*acos(x))
    end
    if n==1
        c(n)=vpa(int(subs(Function,x,cos(x)),x,0,pi),16)/pi;
    else
        c(n)=vpa(int(subs(Function,x,cos(x))*cos((n-1)*x),x,0,pi),16)*2/pi;
    end
end
s_C=dot(c(:),T(1:4).');

%Minimized Interplotation Res.
X_t=solve(T(5)==0);
for n=1:order+1
    f(n)=subs(Function,x,X_t(n));
    for i=1:order+1
        lelement(n,i)=(x-X_t(i))/(X_t(n)-X_t(i));
    end
    l(n)=prod(lelement(n,1:n-1))*prod(lelement(n,n+1:order+1));
end
L_em=dot(f(:),l(:).');

%Substitution
F=zeros(step+1,1);
Legendre=zeros(step+1,1);
Chebyshev=zeros(step+1,1);
Lagrange=zeros(step+1,1);
X=(-1:2/step:1);
for n=1:step+1
    F(n)=double(subs(Function,x,X(n)));
    Legendre(n)=double(subs(s_L,x,X(n)));
    Chebyshev(n)=double(subs(s_C,x,X(n)));
    Lagrange(n)=double(subs(L_em,x,X(n)));
end
errLegendre=Legendre-F;
errChebyshev=Chebyshev-F;
errLagrange=Lagrange-F;
%Cleanup
s_L=vpa(s_L);
s_C=vpa(s_C);
L_em=vpa(L_em);
P=vpa(collect(P));
a=vpa(collect(a));
T=vpa(collect(T));
c=vpa(collect(c));
T=vpa(collect(a));
f=vpa(collect(f));
l=vpa(collect(l));
s_L=vpa(collect(s_L));
s_C=vpa(collect(s_C));
L_em=vpa(collect(L_em));

approx_L=double(vpa(int((Function-s_L)^2,x,-1,1)));
approx_C=double(vpa(int((Function-s_C)^2,x,-1,1)));
approx_Lem=double(vpa(int((Function-L_em)^2,x,-1,1)));
