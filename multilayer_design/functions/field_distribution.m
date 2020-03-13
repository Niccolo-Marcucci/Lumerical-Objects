function [z, nz, P, field] = field_distribution(lambda,theta,d,n,r,t,pol)

if length(theta) ~= 1 || length(r) ~= 1 || length(lambda) ~= 1
    error('Inputs should be scalars')
end

K=2*pi/lambda;
c0 = 299792458;
[d,n,dsub,dair] = prepare_multilayer(d,n);

% determining an optimal resolution, 20 points on the thinnest layer
step = min(d(d~=0))/20;
sz = round(sum(d)/step);          % size of z vector
z  = linspace(0,sum(d),sz);
nz = ones(1,sz);

zv =[0; cumsum(d)];
step = z(2)-z(1);
for i=1:sz
    for j=length(zv):-1:2
        if abs(z(i)-zv(j)) <= step/2
            z(i)=zv(j);
            nz(i)=n(j-1);
        elseif (z(i) < zv(j)) && (z(i) >= zv(j-1))
            nz(i)=n(j-1);
        end
     end
end

% Enforce firs and last refractive indeces
nz(1)=n(1);
nz(end)=n(end);

z_sub = -dsub : step : -step;
sz_sub = length(z_sub);
nz_sub = n(1)*ones(1,sz_sub);

z_air = step : step : dair;
sz_air = length(z_air);
nz_air = n(end)*ones(1,sz_air);

theta=theta/180*pi;
beta = n(1)*sin(theta);
costheta_z = sqrt(nz.^2-beta^2)./nz;  

% There are two methods of proceiding:
% - the first one is the TMM applied from left to right, where
%   Eout = T * Ein
%   i.e. field on te right of the interface is equal to matrix T times
%   field on the left (assuming input field comes from the left
%   All vectors and matrices that refer to this metod have suffix t
% - the second one is the TMM applied from right to left, where
%   Ein = D * Eout
%   All vectors and matrices that refer to this method have suffix d
%
% They are absolutely equivalent, but the matrices T and D use
% diffrent fresnel coefficients and the results might differ a little
% bit when propagation is considered on long distances.
%
%
% The code is set to work only on using the T matrix

Et=zeros(2,sz);
Et(1,1)=1;
Et(2,1)=r;
Mt = eye(2);

Et_air=zeros(2,sz_air);
Et_sub=zeros(2,sz_sub);

for i=1:sz-1
    kt=K*nz(i);
    Pt = prop(kt,z(i+1)-z(i),costheta_z(i));
    Mt = Pt*Mt;
    if nz(i) ~= nz(i+1)
        T = Tij(nz(i),nz(i+1),beta,pol);
        Mt = T*Mt;
    end
    
    Et(:,i+1)  = Mt*[1;r];
end  

% Ed=zeros(2,sz);
% Ed(1,end)=t;
% Ed(2,end)=0;
% Md = eye(2);    
% 
% Ed_air=zeros(2,sz_air);
% Ed_sub=zeros(2,sz_sub);
% 
% for i=sz:-1:2
%     if nz(i-1)~=nz(i)
%         D = Dij(nz(i-1),nz(i),beta,pol);
%         Md= D*Md;
%     end
%     kd=K*nz(i-1);
%     Pd = prop(kd,z(i-1)-z(i),costheta_z(i-1)) ;
%     Md = Pd*Md;
%     
%     Ed(:,i-1) = Md*[t;0];
% end

z  = [z_sub,   z, z_air+z(end)];
nz = [nz_sub, nz, nz_air];

Et_air(1,:) = Et(1,end)*exp(+1i*K*nz(end)*costheta_z(end)*z_air);
Et_sub(1,:) = Et(1,1)*exp(+1i*K*nz(1)*costheta_z(1)*z_sub);
Et_sub(2,:) = Et(2,1)*exp(-1i*K*nz(1)*costheta_z(1)*z_sub);
Et=[Et_sub, Et, Et_air];

% Ed_air(1,:) = Ed(1,end)*exp(+1i*K*nz(end)*costheta_z(end)*z_air);
% Ed_sub(1,:) = Ed(1,1)*exp(+1i*K*nz(1)*costheta_z(1)*z_sub);
% Ed_sub(2,:) = Ed(2,1)*exp(-1i*K*nz(1)*costheta_z(1)*z_sub);
% Ed=[Ed_sub, Ed, Ed_air];


% figure
% plot(z,abs(Ed(1,:)+Ed(2,:)).^2,z,abs(Et(1,:)+Et(2,:)).^2,...
%      z,real(nz)*400);
% nicePlot
% ylim([0 1400])

if nargout > 3
    H_correction=nz/c0;

    field.Er = Et(1,:);
    field.El = Et(2,:);

    field.Hr = Et(1,:)*H_correction;
    field.Hl = Et(2,:)*H_correction;
end

if pol == 'p'
    P = abs((Et(1,:)+Et(2,:)).*nz).^2;
else
    P = abs(Et(1,:)+Et(2,:)).^2;
end

end