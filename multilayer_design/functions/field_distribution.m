function [z, nz, P, field] = field_distribution(lambda,theta,d,n,r,t,pol)

K=2*pi/lambda;
c0 = 299792458;
Z0 = 376.7303136;
[d,n,dsub,dair] = prepare_multilayer(d,n);
% d(1) and d(end) are set to zero

% determining an optimal resolution, 20 points on the thinnest layer
step = min(d(d~=0))/20;
sz = round(sum(d)/step);          % size of z vector
z  = linspace(0,sum(d),sz);
nz = ones(1,sz);

step = z(2)-z(1);

% create nz vector
j=1;
zv=cumsum(d);
za=zv(j);
zb=zv(j+1);
for i=2:sz    
    if abs(z(i)-zb) <= step/2
        j=j+1;
        z(i)=zb;
        za=zv(j);
        zb=zv(j+1);
        nz(i)=n(j+1);
    elseif (z(i) > za) && (z(i) < zb)
        nz(i)=n(j+1);
    end
end
% enforce first and last
nz(1)=n(1);
nz(end)=n(end);

% substrate
z_sub = -dsub : step : -step;
sz_sub = length(z_sub);
nz_sub = n(1)*ones(1,sz_sub);

% external medium
z_air = step : step : dair;
sz_air = length(z_air);
nz_air = n(end)*ones(1,sz_air);

if nargout < 3
    z  = [z_sub,   z, z_air+z(end)];
    nz = [nz_sub, nz, nz_air];
    return
elseif length(theta) ~= 1 || length(r) ~= 1 || length(lambda) ~= 1
    error('Inputs should be scalars')
end

theta=theta/180*pi;
beta = n(1)*sin(theta);
costheta = sqrt(n.^2-beta^2)./n;

% TMM applied from left to right, where
%   Eout = T * Ein
%   i.e. field on te right of the interface is equal to matrix T times
%   field on the left (assuming input field comes from the left
% 

E = zeros(2,sz);
E_air = zeros(2,sz_air);
E_sub = zeros(2,sz_sub);

j=1;
zv=cumsum(d);
za=zv(j);
zb=zv(j+1);

E(:,1) = [1;r];
Mt = Tij(n(1),n(2),beta,pol);

for i=2:sz
    kz = K*n(j+1)*costheta(j+1); 
    Pt = [exp(+1i*kz*(z(i)-za)), 0   ;  
          0 , exp(-1i*kz*(z(i)-za)) ];
    if z(i) == zb
        j=j+1;
        T = Tij(n(j),n(j+1),beta,pol);
        Mt = T*Pt*Mt;
        za = zv(j);
        zb = zv(j+1);
        E(:,i)  = Mt*[1;r];  
    else
        E(:,i)  = Pt*Mt*[1;r];    
    end
end

z  = [z_sub,   z, z_air+z(end)];
nz = [nz_sub, nz, nz_air];

E_air(1,:) = E(1,end)*exp(+1i*K*nz(end)*costheta(end)*z_air);
E_sub(1,:) = E(1,1)*exp(+1i*K*nz(1)*costheta(1)*z_sub);
E_sub(2,:) = E(2,1)*exp(-1i*K*nz(1)*costheta(1)*z_sub);
E=[E_sub, E, E_air];

if nargout > 3    
    H_correction=nz/Z0;
    
    field.Er = E(1,:);
    field.El = E(2,:);
    
    field.Hr = E(1,:)*H_correction;
    field.Hl = E(2,:)*H_correction;
end

if pol == 'p'
    P = abs(sum(E).*nz).^2/2;
else
    P = abs(sum(E)).^2/2;
end

end
