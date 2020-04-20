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

% This function computes the reflectivity and trasmissivity of a
% dielectric multilayer stack. The multilayer vector has to include
% the substrate (as first element) and the external medium (as lat
% layer). 
% The thickness of these two layers will not matter, since the
% computation will start at the first interface and will end at the
% last one. Indeed the values of the field r and t will be computed at
% such interfaces. 

function [R,r,t] = reflectivity_par (lambda,theta_in,d,n,pol)

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

% transverse wavevector.
beta = n(1)*sin(theta_in); 
parfor i=1:size_T
    costheta_z = sqrt(n.^2-beta(i)^2)./n;
    T11=1;
    T12=0;
    T21=0;
    T22=1;
    for k=1:N_layers-1
        n_i = n(k);
        n_j = n(k+1);
        costheta_i = costheta_z(k);
        costheta_j = costheta_z(k+1);
        kz = K*n_i*costheta_i; 
        Pr = exp(+1i*kz*d(k));  
        Pl = exp(-1i*kz*d(k));
        % Fresnel coefficients
        if pol == 's'
            rji = (n_j*costheta_j-n_i*costheta_i)/...
                  (n_i*costheta_i+n_j*costheta_j);
            tji = rji + 1;
        elseif pol == 'p'
            rji = (n_i*costheta_j-n_j*costheta_i)/...
                  (n_j*costheta_i+n_i*costheta_j);
            tji = (rji + 1)*n_j/n_i;
        else
           error("Invalid Polarization. Valid options are 's' or 'p'")
        end

        % expicit matrix product for speed.
        rtij = rji/tji;
        T11t = Pr/tji*T11 + rtij*Pl*T21;
        T12t = Pr/tji*T12 + rtij*Pl*T22;
        T21t = rtij*Pr*T11 + Pl/tji*T21;
        T22t = rtij*Pr*T12 + Pl/tji*T22;
        T11 = T11t;
        T12 = T12t;
        T21 = T21t;
        T22 = T22t;
        
        % next is kept as a reminder
        % T=Tij*P*T; 
        % D=D*P*Dijc;
    end
    r(i) = -T21/T22;
    t(i) = T11+r(i)*T12;
end
R=abs(r).^2;
end
