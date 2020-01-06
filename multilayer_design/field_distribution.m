function [P, z] = field_distribution(lambda,theta,d,n,r,pol)


sizeZ=5e3;
z=linspace(d(1),sum(d),sizeZ);
nz = n(1)*ones(1,sizeZ);
for i=1:sizeZ
    for j=2:length(d)
        if (z(i)<=sum(d(1:j))&&(z(i)>sum(d(1:j-1))))
            nz(i)=n(j);
        end
    end
end

E=zeros(2,sizeZ);

theta_z=asin(n(1)./nz.*sin(theta));

E(1,1)=1;
E(2,1)=r;
for i=1:sizeZ-1
    k=2*pi/lambda*nz(i);
    E(:,i+1)=prop(k,z(i+1)-z(i),theta_z(i))*E(:,i);
    if nz(i)~=nz(i+1)
        E(:,i+1)=Tij(nz(i),nz(i+1),theta_z(i),theta_z(i+1),pol)*E(:,i+1);
    end
end
if pol == 'p'
    H_correction=nz/physconst('lightspeed');
else
    H_correction=ones(1,sizeZ);
end
P = abs((E(1,:)+E(2,:)).*H_correction).^2;

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
    P = [exp(-1i*(kz*d)), 0   ;  
          0 , exp(+1i*(kz*d)) ];
end

end