clear, %close all
tic
theta_v = linspace(40,56,1e3);                    % angle vector
lambda_v=linspace(460,630,1e3)*1e-9;

n_Ti2O2=2.53+1i*1e-4;
n_AlO=1.65+1i*1e-4;
n_SiO2=1.46;
n_pmma=1.48+1i*1e-4;
n_Ta2O5=2.08+1i*1e-4;

n_in=n_SiO2;
n_out=1;
nL=n_SiO2+1i*1e-4;
nH=n_Ta2O5;
ntail=1;%n_pmma;

tail = 75e-9;
dlastSilica = 127e-9;
dL = 137e-9;
dH = 95e-9;

N=7;
% N high-low couples + one for substrate, one for air, one for the pmma layer, one for last titania
% one for allumina protective layer and one for second last titania.
N_layers=2*N+3;    

n1 = zeros(N_layers,1);
d  = zeros(N_layers,1);


d(1) = 1e-6;                        % substrate
d(2:2:2*N) = dH;
d(3:2:2*N+1) = dL;
d(end-2) = dlastSilica;
d(end-1) = tail;
d(end) = 5e-6;


n1(1) = n_in;
n1(2:2:2*N) = nH;
n1(3:2:2*N+1) = nL;
n1(end-1) = ntail;
n1(end) = n_out;

R=zeros(length(lambda_v),length(theta_v));
parfor i=1:length(lambda_v)
    lambda=lambda_v(i);
    [R(i,:),~]=reflectivity(lambda,theta_v,d,n1,'s');
end

[theta,lambda]=meshgrid(theta_v,lambda_v);
% R_p=squeeze(RT1.Rp);
% beta=sin(theta/180*pi)*n_in*2*pi./lambda;
c=physconst('lightspeed');
figure
% s=surf(beta,2*pi*c./lambda,1-R);
s=surf(theta,lambda,1-R);
s.EdgeColor='none';
colormap('hot')
colorbar
view(2)
toc