function [R,r,t,T] = reflectivity (lambda,theta_in,d,n,pol)
theta_in=theta_in/180*pi;
N_layers=length(d);
if N_layers~=length(n) 
    error("thicknesses and refractive index vectors must have same length")
end

size_T=length(theta_in);
r=zeros(1,size_T);
t=zeros(1,size_T);
% A=zeros(2,size_T);

%% loop
k=0;
for i=1:size_T
    theta_z=asin(n(1)./n.*sin(theta_in(i)));
%     D=eye(2);
%     T=eye(2);
    T11=1;
    T12=0;
    T21=0;
    T22=1;
    for j=1:N_layers-1
        P = prop(2*pi/lambda*n(j),d(j),theta_z(j));
%         D=D*P*Dij(n(j),n(j+1),theta_z(j),theta_z(j+1),pol);
        Tijc=Tij(n(j),n(j+1),theta_z(j),theta_z(j+1),pol);
%         T=Tijc*P*T;
        [T11, T12, T21, T22] = product(Tijc*P,T11, T12, T21, T22);
    end
    P = prop(2*pi/lambda*n(end),d(end),theta_z(end));
%     D = D*P;
%     t(i)=1/D(1,1);
%     r(i) = D(2,1)*t(i);
%     T = P*T
%     r(i) = -T(2,1)/T(2,2);
%     A=product(T,[1;r(i)]);
    
    [T11, T12, T21, T22] = product(P,T11, T12, T21, T22);
    r(i) = -cdiv(T21,T22);
    A=product2(T11, T12, T21, T22,[1;r(i)]);
    t(i) = A(1);%*sqrt(n(end)/n(1)*real(cos(theta_z(end)))/cos(theta_z(1)));
    
    if abs(A(2))>1e-13%
        abs(A(2));
        k=k+1; 
    end
end

T=abs(t).^2;
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
    P = [exp(+1i*(kz*d)), 0   ;  
          0 , exp(-1i*(kz*d)) ];
end

function w= product2(M11,M12,M21,M22,v)
    w(1) = cpro(M11,v(1))+cpro(M12,v(2));
    w(2) = cpro(M21,v(1))+cpro(M22,v(2));
end

function [W11, W12, W21, W22]= product(M,V11,V12,V21,V22)
    W11 = cpro(M(1,1),V11)+cpro(M(1,2),V21);
    W12 = cpro(M(1,1),V12)+cpro(M(1,2),V22);
    W21 = cpro(M(2,1),V11)+cpro(M(2,2),V21);
    W22 = cpro(M(2,1),V12)+cpro(M(2,2),V22);
end
function C = cpro(a,b)
    C = abs(a)*abs(b)*exp(1i*(angle(a)+angle(b)));
end
function C = cdiv(a,b)
    C = abs(a)/abs(b)*exp(1i*(angle(a)-angle(b)));
end
end