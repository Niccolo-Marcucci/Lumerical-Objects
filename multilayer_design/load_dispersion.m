clear;
load dispersion
[theta,lambda]=meshgrid(RT1.theta,RT1.lambda);
R_p=squeeze(RT1.Rp);
beta=sin(theta/180*pi)*1.5*2*pi./lambda;
c=physconst('lightspeed');
figure
s=surf(beta,2*pi*c./lambda,1-R_p);
s.EdgeColor='none';
colormap('hot')
colorbar
view(2)