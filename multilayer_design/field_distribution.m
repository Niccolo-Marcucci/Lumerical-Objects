function [P, z, nz] = field_distribution(lambda,theta,d,n,r,t,pol)

if length(theta) ~= 1 || length(r) ~= 1 || length(lambda) ~= 1
    error('Wavelegth, angle of incidence and thus reflectivity value should be scalars')
end
sz=5e3;
z=linspace(0,sum(d),sz);
nz = n(1)*ones(1,sz);
for i=1:sz
    for j=2:length(d)
        if (z(i)<=sum(d(1:j))&&(z(i)>sum(d(1:j-1))))
            nz(i)=n(j);
        end
    end
end

Er=zeros(2,sz);
El=zeros(2,sz);
theta=theta/180*pi;
theta_z=asin(n(1)./nz*sin(theta));

Er(1,1)=1;
Er(2,1)=r;
El(1,end)=t;
El(2,end)=0;
for i=1:sz-1
    k=2*pi/lambda*nz(i);
    Pr = prop(k,z(i+1)-z(i),theta_z(i));
    Pl = prop(k,z(sz-i)-z(sz-i+1)-0,theta_z(sz-i+1)) ;
    Er(:,i+1)  = Pr*Er(:,i);
    El(:,sz-i) = Pl*El(:,sz-i+1);
    if nz(i)~=nz(i+1)
        T = Tij(nz(i),nz(i+1),theta_z(i),theta_z(i+1),pol);
        Er(:,i+1)  = T*Er(:,i+1);
    end
    if nz(sz-i)~=nz(sz-i+1)
        D = Dij(nz(sz-i),nz(sz-i+1),theta_z(sz-i),theta_z(sz-i+1),pol);
        El(:,sz-i) = D*El(:,sz-i+1);
    end
end
if pol == 'p'
    H_correction=nz/physconst('lightspeed');
else
    H_correction=ones(1,sz);
end
figure
plot(z,abs(El(1,:)+El(2,:)),z,abs(El(1,:)),'+',z,abs(El(2,:)));
nicePlot
figure
plot(z,abs(Er(1,:)+Er(2,:)).^2,z,abs(Er(1,:)),'+',z,abs(Er(2,:)));
nicePlot
ylim([0 10])
P = abs((Er(1,:)+El(2,:)).*H_correction).^2;

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
    P = [exp(+1i*(kz*d)), 0   ;  
          0 , exp(-1i*(kz*d)) ];
end

end