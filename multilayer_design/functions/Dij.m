% 
% 
function D = Dij(n_i,n_j,theta_i,theta_j,pol)
    if pol == "s"
        rij = (n_i*cos(theta_i)-n_j*cos(theta_j))./...
              (n_i*cos(theta_i)+n_j*cos(theta_j));
        tij = rij + 1;
    elseif pol == "p"
        rij = (n_j*cos(theta_i)-n_i*cos(theta_j))./...
              (n_j*cos(theta_i)+n_i*cos(theta_j));
        tij = (rij + 1)*n_i/n_j;
    else 
        error("Invalid Polarization. Valid options are 's' or 'p'")
    end

    D = 1/tij* [ 1 , rij ;
                rij, 1  ];  
    
end