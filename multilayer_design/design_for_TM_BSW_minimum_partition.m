theta = linspace(40,60,1000);
lambda=570e-9;

n_AlO=1.65+1i*1e-4;
n_in=1.5;
n_out=1;
nL=1.49+1i*1e-4;
nH=2.53+1i*1e-4;
ntail=1.49;

N=11;

kL=2*pi/lambda*real(nL)*cos(40/180*pi);
kH=2*pi/lambda*real(nH)*cos(40/180*pi);
FF=0.7;
period=pi*0.95/(kL*FF+kH*(1-FF));
dL=period*FF;
dH=period*(1-FF);

dlast=dH*1.45;
dsecondlast=0.2*dH;
dAlO=15e-9;
tail=0;

% N high-low couples + one for substrate, one for air, one for the pmma layer, one for last titania
% one for allumina protective layer and one for second last titania.
N_layers=2*N+6;
value=[0,60e-9];
for l=1:2
    
    tail=value(l);
    
    n1 = ones(1,N_layers);
    d  = zeros(1,N_layers);
    
    n1(1) = n_in;
    n1(2:2:2*N) = nH;
    n1(3:2:2*N+1) = nL;
    n1(end-4) = nH;
    n1(end-3) = n_AlO;
    n1(end-2) = nH;
    n1(end-1) = ntail;
    n1(end) = n_out;
    
    d(1) = 1e-6;
    d(2:2:2*N) = dH;
    d(3:2:2*N+1) = dL;
    d(end-5) = dL*0.5;
    d(end-4) = dsecondlast;
    d(end-3) = dAlO;
    d(end-2) = dlast;
    d(end-1) = tail;
    d(end) = 1e-6;
    
    theta = linspace(40,60,1000);
    [R,r] = reflectivity(lambda,theta,d,n1,'p');
    [~,idx]=findpeaks(1-R,'Npeaks',1,'SortStr','descend');
    
    figure(1)
    plot(theta,R);
    hold
    
    
    [P, z] = field_distribution(lambda,theta(idx),d,n1,r(idx),'p');
    figure(2)
    semilogy(z,P);
end
    
%     %%
%     RT1 = stackrt(n1,d,c/lambda,theta);
%     
%     idx=findpeaks(1-RT1.Rp);
%     field1 = stackfield(n1,d,c/lambda,theta(idx),2e3,0,5e-6);
%     
%     z=field1.z;
%     nz = n1(1)*ones(length(z));
%     for i=1:length(z)
%         for j=2:length(d)
%             if (z(i)<=sum(d(1:j))&&(z(i)>sum(d(1:j-1))))
%                 nz(i)=n1(j);
%             end
%         end
%     end
%     
%     H=abs(pinch(sqrt(field1.Hp(1,1,:,1,1,1)^2+field1.Hp(1,1,:,1,1,2)^2+field1.Hp(1,1,:,1,1,3)^2)))^2;
%     
%     if almostequal(value,0)
%         R1=RT1.Rp;
%         z1=z;
%         H1=H;
%         n_eff1=sin(theta(idx)*pi/180)*n_in;
%     else
%         R2=RT1.Rp;
%         z2=z;
%         H2=H;
%         n_eff2=sin(theta(idx)*pi/180)*n_in;
%     end
%     
% end
% plot(theta,R1,R2,'no PMMA','with PMMA');
% plot(z1,H1,H2,real(nz*300));
% 
% contrast_diff=abs(n_eff2-n_eff1);
% contras_ratio=n_eff2/n_eff1;