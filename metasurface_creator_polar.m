% Copyright 2020 Niccolò Marcucci <niccolo.marcucci@polito.it>
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This code implements the same algorithm proposed in the file
%%% "metasurface_creator.m" code, trying to perform a pseudo extension to
%%% radial symmetries, thus performing the computations in polar
%%% coordinates and then convert to cartesian coordinates.
%%% The code works, but the results were poorly satisfactory, therefore it
%%% contains some incompleted parts and it not weel commented yet. (Yet is
%%% totally similar to the previous code)

clear
close all

eps_av=2;
lambda=0.4e-6;
Rmax = 20e-6;
Rmin = 10e-6;
thetaMax = 4*pi;
K=2*pi/lambda;
Ntheta = ceil((thetaMax*Rmax/lambda)/4);
Nr = ceil((Rmax-Rmin)/lambda*4);
dtheta=thetaMax/Ntheta*Rmin;
dr=(Rmax-Rmin)/Nr;

N_loT=30;
N_lor=ceil(N_loT*9/16);%Nr/Ntheta);
N_hiT=6e3;
N_hiR=ceil(N_hiT*Nr/Ntheta);

theta=linspace(-1,1,Ntheta)*thetaMax/2*Rmin;
r=linspace(Rmin,Rmax,Nr);
[THETA,R]=meshgrid(theta,r);
chirp=R/Rmin;
%[ones(floor(Nr/2),Ntheta);ones(ceil(Nr/2),Ntheta)*Rmax/Rmin]+0.3;

vR=linspace(-0,1,Ntheta);
% b=linspace(-0,1,Nr);
vI=linspace(-0,1,ceil(Nr*0.4));
vI=[zeros(1,floor(floor(Nr*0.6)/2)) vI ones(1,ceil(floor(Nr*0.6)/2))];
[vR,vI]=meshgrid(vR,vI);
gamma=(vI)*pi/2;%angle(a+1i*b);
% gamma(gamma>+pi/2)=gamma(gamma>+pi/2)-pi;
% gamma(gamma<-pi/2)=gamma(gamma<-pi/2)+pi;
% smoother=1./(1+exp(-1.5-3*a));%2*exp(-10*abs(a+1i*b).^2);
% smoother(smoother>1)=1;
% gamma=gamma.*smoother;
% gamma=abs(a+1i*b)*pi/2;
% gamma=cumsum(ones(floor(Nx/2),Ny))'/floor(Nx/2)*pi/2;
% gamma=[zeros(Ny,floor(ceil(Nx/2)/2)) gamma pi/2*ones(Ny,..
%                                           ceil(ceil(Nx/2)/2))];
surface(gamma)
colorbar
colormap('jet');
close all

% write the diff equation grad(phi(x,y))=K(x,y) in the form D *phi = K 
% neuman boundary conditions
N2=Ntheta*Nr;
vec=ones(N2-Nr,1);
vec(1:Nr)=2;
vec=[(1:Nr)'; vec];                     % first N element will be discarded
Dtheta=spdiags([-vec(end:-1:1) vec],[-Nr +Nr],N2,N2);
vec=zeros(N2,1);
vec(1:Nr)=-2;
vec(end-Nr+1:end)=+2;
Dtheta=spdiags(vec,0,Dtheta);


vec=ones(N2-1,1);
vec(Nr:Nr:end)=0;
vec(1:Nr:end)=2;
vec=[0; vec];                          % first N element will be discarded

Dr=spdiags([-vec(end:-1:1) vec],[-1 +1],N2,N2);
vec=zeros(N2,1);
vec(1:Nr:end)=-2;
vec(Nr:Nr:end)=+2;
Dr=spdiags(vec,0,Dr);

Dtheta=Dtheta/(2*dtheta);
Dr=Dr/(2*dr);
% full(Dx)
% full(Dy)

A=([Dtheta;Dr]);
A_=A.'*A;
PHI=zeros(N_hiR,N_hiT,2);
f=figure('units','normalized','outerposition',[0 0 1 1]);

alpha=0;
beta_v=[alpha pi/2+alpha];
for i=1:length(beta_v)
    
    beta=beta_v(i);%(1:10)*pi/20
    ktheta=cos(beta)*K*ones(Nr,Ntheta);
    kr=sin(beta)*K*ones(Nr,Ntheta);
    Ktheta= (cos(gamma).*ktheta - sin(gamma).*kr).*chirp;
    Kr    = (sin(gamma).*ktheta + cos(gamma).*kr);
    
    ktheta=Ktheta(:);
    kr=Kr(:);
    b=([ktheta;kr]);
    b_=A.'*b;
    phi=A_\b_;
    
% %     phi=Dr\kr;
    phi=reshape((phi),Nr,Ntheta);
%     
%     x_=THETA-min(min(THETA));
%     y_=R-min(min(R));
%     phi=(cos(x_./max(max(x_))*pi/2)-sin(y_/max(max(y_))*pi/2))*K;
%     imagesc(theta,r,phi);
%     colorbar;
    newt=linspace(-1,1,N_hiT)*thetaMax/2*Rmin;
    newr=linspace(Rmin,Rmax,N_hiR);
    [newT,newR]=meshgrid(newt,newr);
    PHI(:,:,i)=interp2(THETA,R,phi,newT,newR);
    imagesc(newt,newr,eps_av+(mean(cos(PHI(:,:,i)),3)>0))
    set(gca,'YDir','normal')
    colorbar
    % colormap('gray')
    % xlim([newx(1) newx(end)])
    % ylim([newy(1) newy(end)])
    dR=newr(2)-newr(1);
    dTheta=newt(2)-newT(1);
    [Ktheta_real,Kr_real]=gradient(PHI(:,:,i),dTheta,dR);
    newt2=linspace(-1,1,N_loT)*thetaMax/2*Rmin;
    newr2=linspace(Rmin,Rmax,N_lor);
    [newT2,newR2]=meshgrid(newt2,newr2);
    Ktheta=interp2(THETA,R,Ktheta,newT2,newR2);
    Kr=interp2(THETA,R,Kr,newT2,newR2);
    Ktheta_real=interp2(newT,newR,Ktheta_real,newT2,newR2);
    Kr_real=interp2(newT,newR,Kr_real,newT2,newR2);
    hold on
    
    AR=(f.OuterPosition(3)*f.InnerPosition(3))/...     % 16/9 is the screen
        (f.OuterPosition(4)*f.InnerPosition(4))*16/9;  % aspect ratio
    
    C =range(newr2)/range(newt2)*AR;
    Kmax=max(max(sqrt(Ktheta.^2+Kr.^2)));
    dx=newt2(2)-newt2(1);
    dy=newr2(2)-newr2(1);
    sf=max([dx dy])/Kmax*0.7;
    quiver(newT2,newR2,Ktheta*sf,Kr*sf*C,0,'r','linewidth',3);
    quiver(newT2,newR2,Ktheta_real*sf,Kr_real*sf*C,0,'k','linewidth',3);
%     axis('equal')
    drawnow
    pause(0.5)
    hold off

end

%%
figure('units','normalized','outerposition',[0 0 1 1]);
idx1=(1:ceil(N_hiT/2)+1)+floor(N_hiT/4);
newPHI=PHI(:,idx1,:);
beta=newt(idx1)/Rmin;
r=linspace(Rmin,Rmax,N_hiR);
[BETA,R]=meshgrid(beta,r);
X=R.*cos(BETA);
Y=R.*sin(BETA);

% %%
PHI1=newPHI(:,floor(end/2)+1:end,:);
PHI2=newPHI(:,    1:floor(end/2),:);
dPHI=mean(PHI1(:,end,:)-PHI2(:,1,:));
PHIsmooth=[PHI1 PHI2+dPHI];
% s=surface(newT(:,idx1),newR(:,idx1),(PHIsmooth(:,:,1)));
% s.EdgeColor='none';
% colorbar
idx2=round(length(idx1)*1/3-pi/2):round(length(idx1)*2/3);
for i=1:2
for j=1:N_hiR
PHIsmooth(j,idx2,i)=smooth(PHIsmooth(j,idx2,i),200);
end
PHI1=PHIsmooth(:,ceil(end/2)+1:end,i);
PHI2=PHIsmooth(:,    1:ceil(end/2),i);
newPHI(:,:,i)=[PHI1-dPHI(:,:,i) PHI2];
newdPHI=mean(newPHI(:,end,i)-newPHI(:,1,i));
corr=squeeze(newdPHI-round(newdPHI/2/pi)*2*pi);
CORR=linspace(0,corr,length(beta));
[CORR,~]=meshgrid(CORR,r);
newPHI(:,:,i)=newPHI(:,:,i)-CORR;
end

% newdPHI=mean(newPHI(:,end,:)-newPHI(:,1,:));
% corr=squeeze(newdPHI-round(newdPHI/2/pi)*2*pi)
%
% s=surface(newT(:,idx1),newR(:,idx1),(newPHI(:,:,i)));
% s.EdgeColor='none';
% colorbar
eps=0+(sum(cos(newPHI),3)>0.2);

s=surface(X,Y,1-eps);
s.EdgeColor='none';
axis('equal')
xlim([min(min(X)) max(max(X))]);
ylim([min(min(Y)) max(max(Y))]);
colormap('gray')
set(gca, 'Visible', 'off')
set(gcf, 'units','pixels','outerposition',[0 0 2000 2000],'resize','off')
% save metasurface_ring2 X Y eps
% saveas(figure(1),'metasurface_ring','png')
function imagescG(X,Y,C)

Cg=gpuArray(C);
Xg=gpuArray(X*1e6);
Yg=gpuArray(Y*1e6);

imagesc(Xg,Yg,Cg)
end
