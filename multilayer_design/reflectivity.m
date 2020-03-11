function [R,r,t,A] = reflectivity (lambda,theta_in,d,n,pol)
theta_in=theta_in/180*pi;
N_layers=length(d);
if N_layers~=length(n) 
   error("thicknesses and refractive index vectors must have same length")
end

K=2*pi/lambda;

d(end)=0;
size_T=length(theta_in);
r=zeros(1,size_T);
t=zeros(1,size_T);
A=zeros(2,2,size_T);

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
    A(:,:,i) = T;
%     t(i) = t(i)*sqrt(n(end)/n(1)*real(cos(theta_z(end)))...
%                     /cos(theta_z(1)) );
end

R=abs(r).^2;

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