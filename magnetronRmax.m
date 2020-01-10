function [ rmax ] = magnetronRmax(R1, a, tau, plot)
%MAGNETRONRMAX finds trajectories of magnetron effect electrons
%Inputs: R1 = initial radius of electrons in m
%        a  = parameter (m*V)/(e*B^2*ln(r2/r1)) in m^2
%        tau = nondimentional time interval
%        plot = 'plot' if trajectory plot is desired, omitted otherwise
%Outputs: rmax = max radius of calculated trajectory

    if nargin < 3
        error('too few parameters');
    end
    if nargin < 4 || isempty(plot)      %handle plot as optional parameter
        plot = 'default';
    end
    if ~( strcmpi(plot, 'plot') || strcmpi(plot, 'default'))
        error('invalid input parameter');
    end
    
    dxdt = @(t, x) [x(2);               %define system of DEs using:
                a/x(1) + x(1)*x(4);     %x(1)=r, x(2)=r'
                x(4);                   %x(3)=theta, x(4)=theta'
                -x(2)/x(1)];
     
     options = odeset('RelTol', 1e-6);
    [~, x] = ode45(dxdt, [0 tau], [R1; 0; 0; 0], options);   
    
    r= x(:, 1);   %read radius from first column of returned matrix
    rmax = max(r);
    
    if strcmpi(plot, 'plot')   %if plot is selected
        theta = x(:, 3);        %read angle from third column
        polarplot(theta, r);    %plot trajectory in polar coordinates
    end
end