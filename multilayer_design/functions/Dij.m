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
