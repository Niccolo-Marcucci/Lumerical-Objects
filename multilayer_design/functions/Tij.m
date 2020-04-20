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

% TMM where E_out = T * E_in
% r = -T21/T22, t = T11 + T12*r
function T = Tij(n_i,n_j,beta,pol)
    costheta_i = sqrt(n_i.^2-beta^2)./n_i;  
    costheta_j = sqrt(n_j.^2-beta^2)./n_j;  
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

    T = 1/tji* [ 1 , rji ;
                rji,  1  ];  
    
end
