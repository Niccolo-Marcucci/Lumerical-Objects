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

%%% This trasform a two dimensional binary map (image) into a two
%%% dimensional STL file that can be imported into other software. In order
%%% to extrude in 3D use the software Blender or similar.
%%% Credit to the EM Lab at the University of Texas at El Paso for the idea
%%% of using isocaps (Module 4 of 
%%% https://empossible.net/academics/svl-short-course/)

clear

filename="metasurface_ring2";
load(strcat(filename,'.mat');

% tranform the image into a 3D object by stacking a replica of itself on a
% second layer
eps(:,:,2)=eps;         % esp is the image map
X(:,:,2)=X;             % these are the coordinate. The unit is extreamely  
Y(:,:,2)=Y;             % relevant for further importing
Z=zeros(size(X));       
Z(:,:,2)=1;

tic
[F,V] = isocaps(X,Y,Z,eps,0.5,'zmin');  % extract faces from top or bottom 
toc                                     % layer of the Z map

TR=triangulation(F,V);
% triplot(TR)

stlwrite(TR,strcat(filename,'.stl'),'binary')