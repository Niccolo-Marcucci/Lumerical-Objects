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

%%% This code implements the algorithm proposed in [Raymond C. Rumpf and 
%%% Javier Pazos, "Synthesis of spatially variant lattices," Opt. Express 
%%% 20, 15263-15274 (2012) ] for implementing spatially variant lattices. 
%%% Many thanks to the EM Lab at the University of Texas at El Paso for the
%%% onlie courses and material available at their website 
%%% https://empossible.net/academics/svl-short-course/

clear
close all

eps_av=2;                   % average refractive index
lambda=570e-9/1.1;          % period of the grating

% we can assume that a basilar 2D grating will have two grating vectors 
Ka=2*pi/lambda*[1; 0];      % firs grating vector
Kb=2*pi/lambda*[0; 1];      % second grating vector

Lx = 20e-6;                 
Ly = Lx;
Nx =ceil(Lx/lambda*4);
Ny =ceil(Ly/lambda*4);
dx=Lx/Nx;
dy=Ly/Ny;

N_lox=25;                   % low resolution
N_loy=ceil(N_lox*Ny/Nx);
N_hix=2e3;                  % high resolution
N_hiy=ceil(N_hix*Ny/Nx);

% create a mesh
x=linspace(-Lx,Lx,Nx)/2;
y=linspace(-Ly,Ly,Ny)/2;
[X,Y]=meshgrid(x,y);

%% create arbitrary anle map
v=linspace(0,1,Nx);
w=linspace(0,1,Ny);
[v,w]=meshgrid(v,w);
theta1=angle(v+1i*w);
% next elements are used for creating a more complex angle map
r=sqrt(v.^2+w.^2);
lim1=0.5;
lim2=0.9;
mask =pi/2*(((r>lim1)&(r<lim2)).*((r-lim1)/(lim2-lim1)) +(r>lim2));
% ANGLE MAP
THETA=theta1+ mask;         % set to 0 to ignore it (e.g. when testing)
% plot the map
imagesc(x,y,THETA);
set(gca,'YDir','normal') 
colorbar
colormap('jet');
title('angle map')
% set a breakpoint here to visualize angle map
close all

%% Define finite differene method
% Write the diff equation grad(phi(x,y))=K(x,y) in the form D *phi = K. 
% We assume neuman boundary conditions
N2=Nx*Ny;
% first Dx
vec=ones(N2-Ny,1);
vec(1:Ny)=2;
vec=[(1:Ny)'; vec];                    % first Ny element will be discarded
Dx=spdiags([-vec(end:-1:1) vec],[-Ny +Ny],N2,N2);
vec=zeros(N2,1);
vec(1:Ny)=-2;
vec(end-Ny+1:end)=+2;
Dx=spdiags(vec,0,Dx);
% then Dy
vec=ones(N2-1,1);
vec(Ny:Ny:end)=0;
vec(1:Ny:end)=2;
vec=[0; vec];                          % first Ny element will be discarded
Dy=spdiags([-vec(end:-1:1) vec],[-1 +1],N2,N2);
vec=zeros(N2,1);
vec(1:Ny:end)=-2;
vec(Ny:Ny:end)=+2;
Dy=spdiags(vec,0,Dy);

Dx=Dx/(2*dx);
Dy=Dy/(2*dy);
% compose
A=[Dx;Dy];
A_=A.'*A;

PHI=zeros(N_hiy,N_hix,2);

%% compute PHI for each K-vector of the creating
f=figure('units','normalized','outerposition',[0 0 1 1]);

% in the following, alpha is the angle between the vector Ka (or Kb)
% with respect to the x axis and K is its modulus
K_v=[norm(Ka),norm(Kb)];
alpha_v=[atan2(Ka(2),Ka(1)),atan2(Kb(2),Kb(1))];
for i=1:2
    alpha=alpha_v(i);
    K=K_v(i);
    
    % decompose K vector into Kx,Ky for each point in the grid
    kx=cos(alpha)*K*ones(Ny,Nx);
    ky=sin(alpha)*K*ones(Ny,Nx);

    %%% Modify point by point the value of Kx and Ky
    
    % in this case the vectors are rotated according the the arbitrary
    % angle map THETA
    Kx_new=cos(THETA).*kx - sin(THETA).*ky; % 
    Ky_new=sin(THETA).*kx + cos(THETA).*ky;

    % solve the differential equation
    b=[Kx_new(:);Ky_new(:)];    % transform the matrix into a column vector
    b_=A.'*b;
    phi=A_\b_;
    phi=reshape(phi,Ny,Nx);

    % interpolate on a high resolution mesh
    x_hi=linspace(-Lx,Lx,N_hix)/2;
    y_hi=linspace(-Ly,Ly,N_hiy)/2;
    [X_hi,Y_hi]=meshgrid(x_hi,y_hi);
    PHI(:,:,i)=interp2(X,Y,phi,X_hi,Y_hi);
    
    % plot
    imagesc(x_hi,y_hi,eps_av+cos(PHI(:,:,i)))
    set(gca,'YDir','normal') 
    axis('equal')
    colorbar
    colormap('gray')
    
    % if you want to obseve the shape of phi, uncomment the following
%     s=surf(X,Y,phi);
%     s.EdgeColor='none';
%     colorbar
    
    % compute the resulting local K-vector for the computed phase
    dx=x_hi(2)-x_hi(1);
    dy=y_hi(2)-y_hi(1);
    [Kx_real,Ky_real]=gradient(PHI(:,:,i),dx,dy);
    % create low resolution grid for quiver plot
    x_lo=linspace(-Lx,Lx,N_lox)/2;
    y_lo=linspace(-Ly,Ly,N_loy)/2;
    [X_lo,Y_lo]=meshgrid(x_lo,y_lo);
    % interpolate on such grid both the original and the reuslting Kx,Ky
    Kx_real=interp2(X_hi,Y_hi,Kx_real,X_lo,Y_lo);
    Ky_real=interp2(X_hi,Y_hi,Ky_real,X_lo,Y_lo);
    Kx_lo=interp2(X,Y,Kx_new,X_lo,Y_lo);
    Ky_lo=interp2(X,Y,Ky_new,X_lo,Y_lo);
    % arrow plots: in red the imosed K in green the resulting k
    hold on
    quiver(X_lo,Y_lo,Kx_lo,Ky_lo,'r','linewidth',3);
    quiver(X_lo,Y_lo,Kx_real,Ky_real,'g','linewidth',2);
    drawnow
    pause(0.5)
    hold off
end

%% final plot
f2=figure('units','normalized','outerposition',[0 0 1 1]);


imagesc(x_hi,y_hi,eps_av+(0.5*sum(cos(PHI),3)>0))
set(gca,'YDir','normal') 
axis('equal')
colorbar
colormap('gray')
xlabel('x')
xlabel('y')
title("resulting gratig")

return % run the following only as subsection

%% compose quadrants
% the following is meant at creating a bigger mesh in which the found
% solution is repeted (flipped) on other three quadrants
f3=figure('units','normalized','outerposition',[0 0 1 1]);

newPHI=zeros(2*N_hiy,2*N_hix,2);
for i=1:2
    newPHI(:,:,i)=[PHI(end:-1:1,end:-1:1,i) PHI(end:-1:1,:,i);...
                   PHI(:,end:-1:1,i) PHI(:,:,i)];
end
newPHI=gpuArray(newPHI);
x_hi=gpuArray(linspace(-Lx,Lx,2*N_hix));
y_hi=gpuArray(linspace(-Ly,Ly,2*N_hiy));
[X_hi,Y_hi]=meshgrid(x_hi,y_hi);

% s=surf(x_hi,y_hi,newPHI(:,:,2));
% s.EdgeColor='none';

imagesc(x_hi,y_hi,mean(cos(newPHI(:,:,:)),3)>.1)
set(gca,'YDir','normal') 
axis('equal')
colorbar



