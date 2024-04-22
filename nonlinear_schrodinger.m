% Example usage:
S = 4*pi;
n = 9;
ns = 2^n;
T = 1;
nt = 1000   ;



    % Parameters
    ds = S / (ns-1); % Spatial step size
     
    dt = T / (nt-1);   % Temporal step size
    i = 1i;          

    
    k = 1;
    domain = linspace(-S/2, S/2, ns);

    psi0 = exp(-domain.^2 + i*domain);
    

    % Finite difference method
   
        
        
    kvec = [0:ns/2, -ns/2+1:-1]';
            
            
    
    ode_func = @(t, y) i*((1/2)*((abs(y).^2).*y + ifft( -(kvec.^2) .* fft(y))));
    tspan = [0, T];
    
    [t, y] = ode45(ode_func, tspan, psi0);
    
    

    for j = 1: length(t)
        hold off;
        plot(domain, imag(y(j,:)))
        drawnow



    end

            
    
    % Plot the evolution of the solution in the complex plane
   

    
   

