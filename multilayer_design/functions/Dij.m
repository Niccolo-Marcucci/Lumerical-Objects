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

% TMM where E_in = D * E_out
% t = 1/D11, r = D12*t 
function D = Dij(n_i,n_j,beta,pol)
    costheta_i = sqrt(n_i.^2-beta^2)./n_i;  
    costheta_j = sqrt(n_j.^2-beta^2)./n_j;  
    if pol == "s"
        rij = (n_i*costheta_i-n_j*costheta_j)./...
              (n_i*costheta_i+n_j*costheta_j);
        tij = rij + 1;
    elseif pol == "p"
        rij = (n_j*costheta_i-n_i*costheta_j)./...
              (n_j*costheta_i+n_i*costheta_j);
        tij = (rij + 1)*n_i/n_j;
    else 
        error("Invalid Polarization. Valid options are 's' or 'p'")
    end

    D = 1/tij* [ 1 , rij ;
                rij, 1  ];  
    
end
