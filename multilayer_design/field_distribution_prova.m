function [P, z, nz] = field_distribution(lambda,theta,d,n,r,t,pol)

if length(theta) ~= 1 || length(r) ~= 1 || length(lambda) ~= 1
    error('Wavelegth, angle of incidence and thus reflectivity value should be scalars')
end

K=2*pi/lambda;

dair=d(end);
d(end)=0;

sz=5e3;
z=linspace(0,sum(d),sz);
nz = ones(1,sz);
for i=1:sz
    for j=1:length(d)
        
        % the following is useful to include inside the z vector the 
        % point corresponding to the interfaces
%         if abs(z(i)-sum(d(1:j))) <= z(2)-z(1)/2
%             z(i)=sum(d(1:j));
%         end
        if (z(i)<=sum(d(1:j))&&(z(i)>sum(d(1:j-1))))
            nz(i)=n(j);
        end
    end
end
z_air = z(2)-z(1) : z(2)-z(1) : dair;
sz_air = length(z_air);
nz_air = n(end)*ones(1,sz_air);

theta=theta/180*pi;
theta_z=asin(n(1)./nz*sin(theta));

% There are two methods of proceiding:
% - the first one is the TMM applied from left to right, where
%   Eout = T * Ein
%   i.e. field on te right of the interface is equal to matrix T times
%   field on the left (assuming input field comes from the left
%   All vectors and matrices that refer to this metod have suffix t
% - the second one is the TMM applied from right to left, where
%   Ein = D * Ein
%   All vectors and matrices that refer to this method have suffix d
%
% They are absolutely equivalent, but the matrices T and D use
% diffrent fresnel coefficients.



Et=zeros(2,sz);
Et(1,1)=1;
Et(2,1)=r;
Mt = eye(2);
Et_air=zeros(2,sz_air);

Ed=zeros(2,sz);
Ed(1,end)=t;
Ed(2,end)=0;
Md = eye(2);
Ed_air=zeros(2,sz_air);

for i=1:sz-1
    if nz(sz-i)~=nz(sz-i+1)
        D = Dij(nz(sz-i),nz(sz-i+1),theta_z(sz-i),theta_z(sz-i+1),pol);
        Md= D*Md;
    end
    kd=K*nz(sz-i+1);
    Pd = prop(kd,z(sz-i)-z(sz-i+1),theta_z(sz-i+1)) ;
    Md = Pd*Md;
    
    Ed(:,sz-i) = Md*Ed(:,end);
end

for i=1:sz-1
    kt=K*nz(i);
    Pt = prop(kt,z(i+1)-z(i),theta_z(i));
    Mt  = Pt*Mt;
    
    if nz(i)~=nz(i+1)
        j = 0;
        T = Tij(nz(i),nz(i+1),theta_z(i),theta_z(i+1),pol);
        Mt = T*Mt;
    end
    
    Et(:,i+1)  = Mt*Et(:,1);
end  
z  = [z, z_air+sum(d)];
nz = [nz, nz_air];

Et_air(1,:) = Et(1,end)*exp(+1i*K*nz(end)*cos(theta_z(end))*z_air);
Et=[Et, Et_air];

Ed_air(1,:) = Ed(1,end)*exp(+1i*K*nz(end)*cos(theta_z(end))*z_air);
Ed=[Ed, Ed_air];

if pol == 'p'
    H_correction=nz/physconst('lightspeed');
else
    H_correction=ones(1,sz+sz_air);
end

figure
plot(z,abs(Ed(1,:)+Ed(2,:)).^2,z,abs(Et(1,:)+Et(2,:)).^2,...
     z,real(nz)*400);
nicePlot
ylim([0 1400])
P = abs((Ed(1,:)+Ed(2,:)).*H_correction).^2;

%% functions
function D = Dij(n_i,n_j,theta_i,theta_j,pol)
    if pol == "s"
        rij = (n_i*cos(theta_i)-n_j*cos(theta_j))./...
              (n_i*cos(theta_i)+n_j*cos(theta_j));
        tij = rij + 1;
    elseif pol == "p"
        rij = (n_j*cos(theta_i)-n_i*cos(theta_j))./...
              (n_j*cos(theta_i)+n_i*cos(theta_j));
        tij = (rij + 1)*n_i/n_j;
    else 
        error("Invalid Polarization. Valid options are 's' or 'p'")
    end

    D = 1/tij* [ 1 , rij ;
                rij, 1  ];  
    
end
function T = Tij(n_i,n_j,theta_i,theta_j,pol)
    if pol == 's'
        rij = (n_i*cos(theta_i)-n_j*cos(theta_j))./...
              (n_i*cos(theta_i)+n_j*cos(theta_j));
        rji = -rij;
        tji =  rji + 1;
    elseif pol == 'p'
        rij = (n_j*cos(theta_i)-n_i*cos(theta_j))./...
              (n_j*cos(theta_i)+n_i*cos(theta_j));
        rji = -rij;
        tji = (rji + 1)*n_j/n_i;
    else 
        error("Invalid Polarization. Valid options are 's' or 'p'")
    end

    T = 1/tji* [ 1 , rji ;
                rji,  1  ];  
    
end
function P = prop(k,d,theta)
    kz=k*cos(theta); 
    P = [exp(+1i*kz*d), 0   ;  
          0 , exp(-1i*kz*d) ];
end

end