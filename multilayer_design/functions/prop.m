% 
% 
function P = prop(k,d,theta)
    kz=k*cos(theta); 
    P = [exp(+1i*kz*d), 0   ;  
          0 , exp(-1i*kz*d) ];
end