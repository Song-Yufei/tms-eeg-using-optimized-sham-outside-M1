function V=leadfield1(X,Y,Q,rad,sig,nmax)
% In a layered sphere the program returns the potential V, at  points
% given in X and lying on the outer surface of the sphere, due 
% to  current dipoles  located at points in Y with moments in Q;
% points in Y must lie in the innermost  region of  the sphere.
%
% INPUT: X = (3,m) matrix of field points on the outer spherical surface, 
%        Y = (3,n) matrix of  the locations of the dipoles, 
%        Q = (3,n) matrix of the of moments of the dipoles,
%        rad = the vector of the radii of the shperical regions, 
%            rad(1) < rad(2) < ... < rad(M), M=length(r),
%        sig = the vector of the conductivities sig(j) in the j:th
%              region for j=1:M,
%        nmax = the number of terms in the multipole expansion of the
%               potential; if norm(Y) is not very close to r(1), nmax=30 is
%               usually O.K.
%
% OUTPUT: V = (m,n) matrix containing the potentials: V(j,k) = potential
%              at X(:,k) due to dipole (Y(:,j),Q(:,j)), k=1:m, j = 1:n,
% .........................................................................
% 15 August 2013 : Jukka Sarvas, BECS, Aalto university  
% .........................................................................

m=size(X,2);
n=size(Y,2);
Xhat=X./(ones(3,1)*sqrt(sum(X.^2)));
normY=sqrt(sum(Y.^2));
Yhat=Y./(ones(3,1)*normY);

alpha=ones(m,1)*sum(Q.*Yhat);
beta=Xhat'*Q;
gamma=Xhat'*Yhat;
[L,dL]=legen(nmax,gamma(:));
P=L(2:nmax+1,:);
dP=dL(2:nmax+1,:);
gam=potgam(rad,sig,nmax);

a=0;
b=0;
for k=1:nmax
    c=gam(k)*(ones(m,1)*(normY.^(k-1)));
    a=a+k*c.*reshape(P(k,:)',m,n);
    b=b+c.*reshape(dP(k,:)',m,n);
end
V=1/(4*pi*sig(end))*(alpha.*a+(beta-alpha.*gamma).*b);

        
function [L,dL]=legen(N,x);
% The function returns the values of the Legendre polynomial P_n(x) of
% order n at x in L(n+1,:) for n=0:N .
% If nargout=2, also the derivatives P'_n(x) are returned in dL(n+1,:).
% Note the the legendre functions P_n_1(x) are obtained from P'_n(x) as
%            P_n_1(x)=-sqrt(1-x^2)*P'_n(x)

L=ones(N+1,length(x));
x=x(:).';
if N>0
   L(2,:)=x;
end
if N<2
   L=L(1:N+1,:);
else
   for n=1:N-1
      L(n+2,:)=(2*n+1)/(n+1)*x.*L(n+1,:)-n/(n+1)*L(n,:);
   end
end

if nargout==2
    dL=zeros(N+1,length(x));
    dL(2,:)=ones(1,length(x));
    for n=1:N-1
        dL(n+2,:)=dL(n,:)+(2*n+1)*L(n+1,:);
    end
end

function gam=potgam(r,sig,nmax);
% The function returns the translation coefficients gam(n), 
% INPUT: r = vector of (outer) radii of the layers,
%        sig = vector of the conductivities of the layers,
%        nmax = the maximal n index of the multipole expansion,
% OUTPUT: gam = the column vector of the translation coefficients,

M=length(r);
gam=zeros(nmax,1);
s=sig;
for n=1:nmax
    C=eye(2);
    for j=1:M-1
        c=1/((2*n+1)*s(j+1));
        Cj=c*[n*s(j)+(n+1)*s(j+1),(n+1)*(s(j+1)-s(j))/r(j)^(2*n+1);...
              n*r(j)^(2*n+1)*(s(j+1)-s(j)),(n+1)*s(j)+n*s(j+1)];
        C=Cj*C;
    end
    a=(2*n+1)/(n+1)*r(M)^n;
    b=-C(2,1)+n/(n+1)*r(M)^(2*n+1)*C(1,1);
    gam(n)=a/b;
end


