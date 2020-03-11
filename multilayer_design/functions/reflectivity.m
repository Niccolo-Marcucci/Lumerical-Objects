function [R,r,t] = reflectivity (lambda,theta_in,d,n,pol)

theta_in=theta_in/180*pi;
N_layers=length(d);
if N_layers~=length(n) 
   error("thicknesses and refractive index vectors must have same length")
end

K=2*pi/lambda;

i=length(d);
while n(i-1)==n(end)
    d(i-1)=0;
    i=i-1;
end
% d(1)=0;
d(end)=0;
size_T=length(theta_in);
r=zeros(1,size_T);
t=zeros(1,size_T);

%% loop
for i=1:size_T
    theta_z=asin(n(1)./n.*sin(theta_in(i)));
%     D=eye(2);
    T=eye(2);
    for j=1:N_layers-1
        P = prop(K*n(j),d(j),theta_z(j));
%         D=D*P*Dij(n(j),n(j+1),theta_z(j),theta_z(j+1),pol);
        Tijc=Tij(n(j),n(j+1),theta_z(j),theta_z(j+1),pol);
        T=Tijc*P*T;
    end
%     t(i) = 1/D(1,1);
%     r(i) = D(2,1)*t(i);
    r(i) = -T(2,1)/T(2,2);
    t(i) = T(1,1)+r(i)*T(1,2);
%     t(i) = t(i)*sqrt(n(end)/n(1)*real(cos(theta_z(end)))...
%                     /cos(theta_z(1)) );
end

R=abs(r).^2;

end