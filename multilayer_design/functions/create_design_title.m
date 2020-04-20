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

% useful only if we assume the disign is in the form 
% Substrate-(A - B)*N-X-X-X-...-External_Medium

function title = create_design_title(design_file)
load(design_file,'d_layers','idx_layers');
d_layers=round(d_layers*1e9);
N=length(d_layers);
idx_dA=(d_layers==d_layers(2));
idx_dB=(d_layers==d_layers(3));
if ( idx_layers(idx_dA)~=idx_layers(2)) | ...
                (idx_layers(idx_dA)~=idx_layers(2) )
    error("this multilayer cannot be used with this function")
end

Na=sum(idx_dA);
Nb=sum(idx_dB);
if Na~=Nb
    N_couples=min(Na,Nb);
else
    N_couples=Na;
end

title = [ strcat("", idxToStr(idx_layers(1))," -(",...
                 idxToStr(idx_layers(2)),"-",...
                 idxToStr(idx_layers(3)),")x",string(N_couples));...
          strcat("Thicknesses:  ~ -(",...
                 string(d_layers(2)),"-",...
                 string(d_layers(3)),")x",string(N_couples))];
for i=1+2*N_couples+1:N-1
    title = [strcat(title(1),"- ",idxToStr(idx_layers(i)));...
             strcat(title(2),"- ",string(d_layers(i)))];
end

title(2)=strcat(title(2),"  nm");
end

% has to be updated to the expected meterial library
function str = idxToStr(value)
switch real(value)
    case 1.46
        str = "SiO_2";
    case 1.65
        str = "Al_2O_3";
    case 2.53
        str = "TiO_2";
    case 1.48
        str = "PMMA";
    case 2.08
        str = "Ta_2O_5";
    case 1
        str = "Air";
    otherwise
        str = strcat("fake",string(real(value)),"");
end
end

        