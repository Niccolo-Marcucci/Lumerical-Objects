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

% Computes the field distribution inside a dielectric multilayer
% stack. 
% - The output is computed on the normalized H_tilde=Z0*H;
% - The reolution parameter can be neglected
% - The function can be used also to etract only nz. in this case use
%   in the form:
%   [z, nz] = field_distribution(lambda,theta,d,n,'',res)
function [z, nz, P, field] = field_distribution(lambda,theta,d,n,pol,res)

[d,n,dsub,dair] = prepare_multilayer(d,n);  % check fun description 
                                            % for more details.                                      
% determining an optimal resolution, 50 points on shortest wavelegth
if nargin < 6
    res = 50;
end
step = lambda/max(real(n))/res;
sz = round(sum(d)/step);          % size of z vector
z  = linspace(0,sum(d),sz);
nz = ones(1,sz);
step = z(2)-z(1);

% create nz vector: if a point is located exactly at an interface, it
% will be cosidered to belong to the layer on the right hand side of
% the interface
j=1;
zv=cumsum(d);
za=zv(j);
zb=zv(j+1);
for i=2:sz    
    if abs(z(i)-zb) <= step/2
        j=j+1;
        z(i)=zb;
        za=zv(j);
        zb=zv(j+1);
        nz(i)=n(j+1);
    elseif (z(i) > za) && (z(i) < zb)
        nz(i)=n(j+1);
    end
end
% enforce first and last
nz(1)=n(1);
nz(end)=n(end);

% substrate
z_sub = -dsub : step : -step;
sz_sub = length(z_sub);
nz_sub = n(1)*ones(1,sz_sub);

% external medium
z_air = step : step : dair;
sz_air = length(z_air);
nz_air = n(end)*ones(1,sz_air);

if nargout < 3
    z  = [z_sub,   z, z_air+z(end)];
    nz = [nz_sub, nz, nz_air];
    return
elseif length(theta) ~= 1 || length(lambda) ~= 1
    error('Inputs should be scalars')
end 

% TMM applied from left to right, where
%   Eout = T * Ein
%   i.e. field on te right of the interface is equal to matrix T times
%   field on the left (assuming input field comes from the left)
% Check the appendix for more info.

% first of all extract the reflected field at the first interface
[~,r] = reflectivity(lambda,theta,d,n,pol); 

% determine wave direction in each layer
K=2*pi/lambda;
theta = theta/180*pi;
beta = n(1)*sin(theta);
costheta = sqrt(n.^2-beta^2)./n;  

% and now propagate the fields
E = zeros(2,sz);
E_air = zeros(2,sz_air);
E_sub = zeros(2,sz_sub);

j=1;
zv=cumsum(d);
za=zv(j);
zb=zv(j+1);

E(:,1) = [1;r];
Mt = Tij(n(1),n(2),beta,pol);

for i=2:sz
    kz = K*n(j+1)*costheta(j+1); 
    Pt = [exp(+1i*kz*(z(i)-za)), 0   ;  
          0 , exp(-1i*kz*(z(i)-za)) ];
    if z(i) == zb
        j=j+1;
        T = Tij(n(j),n(j+1),beta,pol);
        Mt = T*Pt*Mt;
        za = zv(j);
        zb = zv(j+1);
        E(:,i)  = Mt*[1;r];  
    else
        E(:,i)  = Pt*Mt*[1;r];    
    end
end

z  = [z_sub,   z, z_air+z(end)];
nz = [nz_sub, nz, nz_air];

E_air(1,:) = E(1,end)*exp(+1i*K*nz(end)*costheta(end)*z_air);
E_sub(1,:) = E(1,1)*exp(+1i*K*nz(1)*costheta(1)*z_sub);
E_sub(2,:) = E(2,1)*exp(-1i*K*nz(1)*costheta(1)*z_sub);
E=[E_sub, E, E_air];


% apply Maxwell equations in order to extract H and E
costheta_z = sqrt(nz.^2-beta^2)./nz;
if pol == 'p'
    % In order to determine Ey and Ez it is necessary to take into
    % account the conventions of signs with which the Fresnel 
    % coefficient have been analytically derived. 
    % Check the appendix for more info.
    sintheta_z = beta./nz;
    
    field.Ery = E(1,:).*costheta_z;
    field.Ely = E(2,:).*costheta_z;
    field.Erz =-E(1,:).*sintheta_z;
    field.Elz =+E(2,:).*sintheta_z;
    field.Ey = field.Ery - field.Ely;
    field.Ez = field.Erz - field.Elz;
    
    field.Hrx = beta*field.Erz - nz.*costheta_z.*field.Ery;
    field.Hlx =-beta*field.Elz - nz.*costheta_z.*field.Ely; 
    field.Hx = field.Hrx + field.Hlx;
    
    P = abs(field.Hx).^2/2;
else
    
    % With s-polarisation the sign convention is much easier. We
    % simply apply the third Maxwell equation
    % Check the appendix for more info.
    
    field.Erx = E(1,:);
    field.Elx = E(2,:);
    field.Ex = field.Erx + field.Elx;
    
    field.Hry = nz.*costheta_z.*field.Erx;
    field.Hly =-nz.*costheta_z.*field.Elx;
    field.Hrz =-beta*field.Erx;
    field.Hlz =-beta*field.Elx;
    
    field.Hy = field.Hry + field.Hly;
    field.Hz = field.Hrz + field.Hlz;
    
    P = abs(field.Ex).^2/2;
end
end
