function [r,t,x,Nn,E] = tmm(wvl,fulln, fullw)
% Transfer matrix method
%
% Inputs
% wvl: wavelength
% fulln: array with all refractive indices
% fullw: array with all layer widths
%
% Outputs:
% r: reflection coefficient
% t: transmission coefficient
% x: position
% Nn: refractive indices (vs x)
% E: optical field (vs x)

    % Exponent factors
    d = fullw*2*pi/wvl.*fulln;

    % Initiate arrays
    x = [];
    E = [];
    Nn = [];
    N = length(fulln);

    % Loop through layers
    for n=1:(N-1)

        % n of adjacent layers
        n1 = fulln(n);
        n2 = fulln(n+1);

        % Fresnel relations
        rs(n) = (n1 - n2)/(n1+n2);
        ts(n) = 2*n1/(n1+n2);

        % Compose transfer matrix
        M(:,:,n) =  [exp(-1i*d(n)),0;0,exp(1i*d(n))] * [1, rs(n); rs(n),1] * 1/ts(n);

        % Multiply with full matrix (if exists)
        if n>=2
            Mt = Mt*M(:,:,n);
        else
            Mt = M(:,:,1);
        end

    end

    % Reflection and transmission coefficients
    r = Mt(2,1)/Mt(1,1);
    t = 1/Mt(1,1);

    % Initiate arrays
    v1 = zeros(1,length(fullw));
    v2 = zeros(1,length(fullw));
    v1(1) = 1;
    v2(1) = r;

    % Loop through layers
    for n=2:N

        % Coefficients
        vw = M(:,:,n-1)\[v1(n-1);v2(n-1)];
        v1(n) = vw(1);
        v2(n) = vw(2);

        % Location array
        xloc = 0:5:fullw(n);

        % Electric fields
        Eloc1 = v1(n)*exp(1i*2*pi/wvl*fulln(n)*xloc);
        Eloc2 = v2(n)*exp(-1i*2*pi/wvl*fulln(n)*xloc);

        % Append to arrays
        x = [x,xloc+sum(fullw(1:n-1))];
        E = [E,(Eloc1+Eloc2)];
        Nn = [Nn,fulln(n)+(xloc*0)];

    end

    % Sort arrays
    [x,ix] = unique(x);
    E = E(ix);
    Nn = Nn(ix);

end