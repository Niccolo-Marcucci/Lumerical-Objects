% This function computes the reflectivity and trasmissivity of a
% dielectric multilayer stack. The multilayer vector has to include
% the substrate (as first element) and the external medium (as lat
% layer). 
% The thickness of these two layers will not matter, since the
% computation will start at the first interface and will end at the
% last one. Indeed the values of the field r and t will be computed at
% such interfaces. 



function [R,r,t,Tr] = reflectivity (lambda,theta_in,d,n,pol)

N_out = nargout;
N_layers=length(d);
if N_layers~=length(n) 
   error("thicknesses and refractive index vectors must have same length")
end

% if you propose a multilayer which is altready optimized for this
% calculation (i.e. no dummy layers and no zero thicknes layers), you
% can comment the next line and save some computational time.
% [d,n] = prepare_multilayer(d,n);

N_layers = length(d);

theta_in = theta_in/180*pi;

K = 2*pi/lambda;

size_T = length(theta_in);
r = zeros(1,size_T);
t = zeros(1,size_T);
Tr= zeros(1,size_T);

beta = n(1)*sin(theta_in);
%% loop
parfor i=1:size_T
    costheta_z = sqrt(n.^2-beta(i)^2)./n;
    T11=1;
    T12=1;
    T21=1;
    T22=1;
%     T = eye(2);  
%     D = eye(2);
    for k=1:N_layers-1
%         [Pr,Pl] = prop(K*n(j),d(j),costheta_z(j));
        n_i=n(k);
        n_j=n(k+1);
        costheta_i=costheta_z(k);
        costheta_j=costheta_z(k+1);
        kz = K*n_i*costheta_i; 
        Pr = exp(+1i*kz*d(k));  
        Pl = exp(-1i*kz*d(k));
%         [T11c,T12c,T21c,T22c] = ...
%                  Tijx(n(j),n(j+1),costheta_z(j),costheta_z(j+1),pol);
        if pol == 's'
            rij = (n_i*costheta_i-n_j*costheta_j)./...
                  (n_i*costheta_i+n_j*costheta_j);
            rji = -rij;
            tji =  rji + 1;
        elseif pol == 'p'
            rij = (n_j*costheta_i-n_i*costheta_j)./...
                  (n_j*costheta_i+n_i*costheta_j);
            rji = -rij;
            tji = (rji + 1)*n_j/n_i;
        else
           error("Invalid Polarization. Valid options are 's' or 'p'")

        end

%         T11c = 1/tji;
        T12c = rji/tji;
%         T21c = T12c;
%         T22c = T11c;
        
        T11=Pr*(T11/tji+T12c*T21);
        T12=Pr*(T12/tji+T12c*T22);
        T21=Pl*(T12c*T11+T21/tji);
        T22=Pl*(T12c*T12+T22/tji);
        
%         T=Tijc*P*T;
%         P = prop(K*n(j),-d(j),costheta_z(j));
%         Dijc=Dij(n(j),n(j+1),beta_i,pol);
%         D=D*P*Dijc;
    end
%     t(i) = 1/D(1,1);
%     r(i) = D(2,1)*t(i);
    r(i) = -T21/T22;
    t(i) = T11+r(i)*T12;
    if N_out > 3
        Tr(i) = abs( t(i)*n(end)/n(1)*real(costheta_z(end))...
                        /costheta_z(1) )^2;
    end
end

R=abs(r).^2;
end

% function [T11,T12,T21,T22] = Tijx(n_i,n_j,costheta_i,costheta_j,pol)
% 
%     if pol == 's'
%         rij = (n_i*costheta_i-n_j*costheta_j)./...
%               (n_i*costheta_i+n_j*costheta_j);
%         rji = -rij;
%         tji =  rji + 1;
%     elseif pol == 'p'
%         rij = (n_j*costheta_i-n_i*costheta_j)./...
%               (n_j*costheta_i+n_i*costheta_j);
%         rji = -rij;
%         tji = (rji + 1)*n_j/n_i;
%     else 
%         error("Invalid Polarization. Valid options are 's' or 'p'")
%     end
% 
%     T11 = 1/tji;
%     T12 = rji/tji;
%     T21 = T12;
%     T22 = T11;
% %     T = 1/tji* [ 1 , rji ;
% %                 rji,  1  ];  
%     
% end