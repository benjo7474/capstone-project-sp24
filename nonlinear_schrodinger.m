clear; close all;

ns = 2^5;
T = 1;
k = 6*pi;
domain = linspace(0, 1, ns+1);
domain = domain(1:ns);
% psi0 = exp(-40*(domain-0.5).^2 + 2i*pi*domain); psi0 = psi0/trapz(domain, psi0);
psi0 = exp(1i*k*domain);

%% Plot IC
figure; hold on;
plot(domain, real(psi0));
plot(domain, imag(psi0));
legend('Real Part', 'Imag Part')

%% Spectral differentiation with method of lines
[t, y] = ode23(@(t,y)ode_func(t,y,ns), [0 T], psi0);

%% Finite differences with method of lines
e = ones(ns,1);
ds = domain(2)-domain(1);
D2 = spdiags([e -2*e e], -1:1, ns, ns);
D2(1,ns) = 1; D2(ns,1) = 1;
D2 = D2/ds^2;
tspan = [0, T];
ode_func = @(t,y) 1i*(0.5*(abs(y).^2).*y + D2*y);
[t,y] = ode45(ode_func, tspan, psi0);

%% Plot
figure;
exact = @(t,x) exp(1i*(k*x - (k^2-0.5)*t)); % wave speed is still too fast. why?
for j = 1:3:length(t)
    plot(domain, real(y(j,:)), 'LineWidth', 2);
    hold on;
    plot(domain, imag(y(j,:)), 'LineWidth', 2);
    soln = exact(t(j),domain);
    plot(domain, real(soln), 'LineWidth', 2);
    plot(domain, imag(soln), 'LineWidth', 2);
    legend('Real Part (Numerical)', 'Imag Part (Numerical)', 'Real Part (Exact)', 'Imag Part (Exact)');
    ylim([-5 5])
    title(['$t=',num2str(t(j)),'$'],'Interpreter', 'latex');
    drawnow
    hold off;
    pause(0.05)
end

%% Exact Solution
exact = @(t,x) exp(1i*(k*x + (k^2-0.5)*t));
figure;
for j=1:2:length(t)
    soln = exact(t(j),domain);
    plot(domain, real(soln));
    hold on;
    plot(domain, imag(soln));
    title(['$t=',num2str(t(j)),'$'],'Interpreter', 'latex');
    drawnow;
    hold off;
end


% function rhs = ode_func(t,y,ns)
% 
% % Make sure y is a column vector
% y = y(:);
% kvec = [0:ns/2, -ns/2+1:-1]';
% rhs = 1i*(0.5*(abs(y).^2).*y + ifft(-kvec.^2 .* fft(y)));
% 
% end

