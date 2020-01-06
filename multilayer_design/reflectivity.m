function [R,r] = reflectivity (lambda,theta_in,d,n,pol)
theta_in=theta_in/180*pi;
N_layers=length(d);
if N_layers~=length(n) 
    error("thicknesses and refractiv index vectors must have same length")
end

size_T=length(theta_in);
r=zeros(1,size_T);

%% loop

for i=1:size_T
    theta_z=asin(n(1)./n.*sin(theta_in(i)));
    D=eye(2);
    for j=1:N_layers-1
        P = prop(2*pi/lambda*n(j),d(j),theta_z(j));
        D=D*P*Dij(n(j),n(j+1),theta_z(j),theta_z(j+1),pol);
    end
    P = prop(2*pi/lambda*n(end),d(end),theta_z(end));
    D = D*P;
    t=1/D(1,1);
    r(i) = D(2,1)*t;
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
function P = prop(k,d,theta)
    kz=k*cos(theta); 
    P = [exp(-1i*(kz*d)), 0   ;  
          0 , exp(+1i*(kz*d)) ];
end

end