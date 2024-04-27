clear all;
close all;
flag = "Eight";

aFlag = true;

N = 500;

t = linspace(0,1,N);

% dt = (t(end)-t(1))/(N-1);

% Defining the curves chosen by the user

switch flag

case "Bezier"

Pts = [0 1 1 2;

0 randi(5)-2 randi(5)-3 0;

0 1 2 3];

Curv = Bezier(Pts, t*2);

case "Helix"

x = cos(pi.*2*t);

y = sin(pi.*2*t);

z = t*2;

Curv = [x; y; z];

case "Viviani Curve"

x = cos(pi.*2.*t).^2;

y = cos(pi.*2.*t).*sin(pi.*2.*t);

z = sin(pi.*2.*t);

Curv = [x; y; z];

case "Torus"

x = 1/3*(cos(pi.*2.*t) + 2*cos(2*pi.*2.*t));

y = 1/3*(sin(pi.*2.*t) - 2*sin(2*pi.*2.*t));

z = 1/3*(2*sin(pi.*2.*t));

Curv = [x; y; z];

case "Kink"

x = 1/3*(cos(pi.*2.*t) + 2*cos(2*pi.*2.*t));

y = 1/3*(sin(pi.*2.*t) - 2*sin(2*pi.*2.*t));

z = 1/3*(2*sin(pi.*2.*t));

Curv = [x; y; z];

case "New Knot"

x = cos(pi*2.*2.*t).*(3+cos(pi*3.*2.*t));

y = sin(pi*3.*2.*t).*(3+cos(pi*3.*2.*t));

z = sin(pi*4.*2.*t);

Curv = [x; y; z];

case "Circle"

x = cos(pi*2.*t);

y = sin(pi*2.*t);

z = zeros(size(t));

Curv = [x; y; z];

case "Circle2"

t=(t.^2+t).*0.5;

x = cos(pi*2.*t);

y = sin(pi*2.*t);

z = zeros(size(t));

Curv = [x; y; z];

case "My Knot"

    t=t.*2.*pi;
        A = 1;          % Amplitude
mu = pi;         % Mean
sigma = 0.2;      % Standard deviation


% Compute Gaussian function
gauss = A * exp(-(t - mu).^2 / (2 * sigma^2));


x=cos(t)+0.5.*gauss;
y=gauss.*-3.*sin(t)+sin(t);
z=0.3.*gauss.*-3.*sin(t);
Curv = [x;y;z];


    case "Cash"

        % Define the helix and motion parameters
radius = 10;            % Radius of the main helix
small_radius = 0.2;     % Radius of the high-frequency corkscrew
pitch = 15;             % Vertical distance for a single complete loop of the helix
high_freq = 15;         % Frequency multiplier for the high-frequency corkscrew
% num_points = 10000;     % Number of points along the helix
% dt = 0.05;              % Time step
% total_time = 10.0;      % Total time for the simulation
% time_steps = total_time / dt;

% Parametrize the helix 
t = t*2*pi;

x = radius * cos(t) - cos(t) .* small_radius .* cos(high_freq * t);
y = radius * sin(t) - sin(t) .* small_radius .* cos(high_freq * t);
z = small_radius * sin(high_freq * t);
Curv = [x;y;z];


    case "Eight"
        a=2*pi*t;
        x = (2 + cos(2*a)).*cos(3*a);
y = (2 + cos(2*a)).*sin(3*a);
z = sin(4*a);
Curv = [x;y;z];

    case "Trefoil"
        a=2*pi*t;
        x = sin(a) + 2*sin(2*a);
y = sin(3*a);
z = cos(a) - 2*cos(2*a);
Curv = [x;y;z];

case "Granny"
        a=2*pi*t;
       x = -0.22*cos(a) - 1.28*sin(a) - 0.44*cos(3*a) - 0.78*sin(3*a);
y = -0.1*cos(2*a) - 0.27*sin(2*a) + 0.38*cos(4*a) + 0.46*sin(4*a);
z = 0.7*cos(3*a) - 0.4*sin(3*a);
Curv = [x;y;z];


end

plot3(x,y,z)

%% Get uniform Spacing

%make uniformly spaced

pt = interparc(300,x,y,z);
pt = pt';
% pt=interparc(100,pt(1,:),pt(2,:),pt(3,:));
% pt=pt';

figure;
hold on
scatter3(pt(1,:),pt(2,:),pt(3,:))
plot3(pt(1,:),pt(2,:),pt(3,:))

edgeVals = edges(pt);

norms = sqrt(sum(edgeVals.^2, 1)); % Calculate norm for each vector

% Plot the norms
figure
plot(norms, 'LineWidth', 2);
xlabel('Index');
ylabel('Norm');
title('Norm of Deriv. Vectors along Curve');
grid on;




%% Compute curvature and torsion

[T,N,B,k,tor] = frenet(pt(1,:),pt(2,:),pt(3,:));
figure

hold on

plot(k, 'DisplayName', 'Curvature (k)');
plot(tor, 'DisplayName', 'Torsion (tor)');

legend;  % Add legend to identify the curves
xlabel('X-axis');  % Label the x-axis
ylabel('Y-axis');  % Label the y-axis
title('Curvature (k) and Torsion (tor)');  % Add a title


%% Compute PSI for NLSE

totalTor = cumsum(tor);

figure

hold on

plot(k, 'DisplayName', 'Curvature (norm of psi)');
plot(tor, 'DisplayName', 'Torsion (tor)');
plot(totalTor, 'DisplayName', 'Total Torsion (argument of psi)');

legend;  % Add legend to identify the curves
xlabel('X-axis');  % Label the x-axis
ylabel('Y-axis');  % Label the y-axis
title('Total Torsion (integral of tor)');  % Add a title

%% TRY ODE45

% Initial curve and time span
% Curv = initial_curve; % Replace initial_curve with your initial curve
y0 = Curv(:); % Reshape the curve to a column vector
tspan = [0,5]; % Time span for simulation

% Solve the ODE using ode45
[t, y] = ode23(@curveODE, tspan, y0);

%% PLOT ANIMATION
% Create a figure


figure
plot(t)
pause(1);
figure;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Curve Animation');

% Initialize animated line
curve_line = animatedline('Color', 'b', 'LineWidth', 2);

num_points = numel(y0) / 3;

view(45, 20);
% Set the axis limits (optionally remove this so it just fits the curves)
lim=3.5;
grid on;
xlim([-lim, lim]);
ylim([-lim, lim]);
zlim([-lim, 10]);

WINDOW=5; %%% WINDOW will plot every window'th iteration, use this to speed up or slow down visualization
for i = 1:WINDOW:numel(t)
    % Reshape the solution to get the points at the current time step
    points_at_current_time = reshape(y(i,:)', [3, num_points]);

    % Clear the previous points
    clearpoints(curve_line);


    % Add the current line
    addpoints(curve_line, points_at_current_time(1,:), points_at_current_time(2,:), points_at_current_time(3,:));

    drawnow;

    % Pause for animation speed
    % pause(0.1);
end



%% PLOT THE CURVE AND BINORMAL VECTORS

[Deriv] = derivFunc(Curv);

hold on;

% Plot the vectors
quiver3(Curv(1,:), Curv(2,:), Curv(3,:), Deriv(1,:), Deriv(2,:), Deriv(3,:), 'r', 'LineWidth', 2);

plot3(Curv(1,:), Curv(2,:), Curv(3,:), 'b'); 

hold off;

% Add labels, title, etc. if needed
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Vector field along curve');
grid on;


%% FUnctions

function [kB]=derivFunc(Curve)

Curve = Curve(:, 1:end-1); %remove last point so overlaps don't fuck things up
    
% Get the number of points
num_points = size(Curve, 2);

% Initialize the edge and derivative matrix
edges = zeros(3, num_points);
kB= zeros(3, num_points);

% Compute edges
for i = 1:num_points
    % Current point
    current_point = Curve(:, i);
    
    % Next point (if we are at the last point, loop back to the first point)
    next_point = Curve(:, mod(i, num_points)+1);
    
    % Compute edge vector
    edge_vector = next_point - current_point;
    
    % Store edge vector in the edges matrix
    edges(:, i) = edge_vector;
end

for i = 1:num_points
    kB(:,i)=2.*cross(edges(:,mod(i-2, num_points)+1),edges(:,i))./(norm(edges(:,mod(i-2, num_points)+1)).*norm(edges(:,i))+dot(edges(:,mod(i-2, num_points)+1),edges(:,i)));%from paper, except further divided by the average norm of the vectors?
    kB(:,i)=kB(:,i)./(norm(edges(:,mod(i-2, num_points)+1))+norm(edges(:,i))).*2;
end

edges = [edges, edges(:, 1)]; %Add on last points back
kB = [kB, kB(:, 1)]; %
% kB=edges; %For testing REMOVE THIS LATER
end

function [edges]=edges(Curve)

Curve = Curve(:, 1:end-1); %remove last point so overlaps don't fuck things up
    
% Get the number of points
num_points = size(Curve, 2);

% Initialize the edge and derivative matrix
edges = zeros(3, num_points);
kB= zeros(3, num_points);

% Compute edges
for i = 1:num_points
    % Current point
    current_point = Curve(:, i);
    
    % Next point (if we are at the last point, loop back to the first point)
    next_point = Curve(:, mod(i, num_points)+1);
    
    % Compute edge vector
    edge_vector = next_point - current_point;
    
    % Store edge vector in the edges matrix
    edges(:, i) = edge_vector;
end

for i = 1:num_points
    kB(:,i)=2.*cross(edges(:,mod(i-2, num_points)+1),edges(:,i))./(norm(edges(:,mod(i-2, num_points)+1)).*norm(edges(:,i))+dot(edges(:,mod(i-2, num_points)+1),edges(:,i)));%from paper, except further divided by the average norm of the vectors?
    kB(:,i)=kB(:,i)./(norm(edges(:,mod(i-2, num_points)+1))+norm(edges(:,i))).*2;
end

edges = [edges, edges(:, 1)]; %Add on last points back
kB = [kB, kB(:, 1)]; %
% kB=edges; %For testing REMOVE THIS LATER
end


% ODE function
function dydt = curveODE(~, y)
    
    % Reshape y back to a matrix of points
    num_points = numel(y) / 3;
    points = reshape(y, [3, num_points]);

    % Compute the derivative of each point using derivFunc
    dydt = derivFunc(points);
    
    % Reshape dydt to a column vector
    dydt = dydt(:);
end





%%
function [T,N,B,k,tor] = frenet(x,y,z)

    if nargin == 2

        z = zeros(size(x));

    end

    % CONVERT TO COLUMN VECTOR

    x = x(:);

    y = y(:);

    z = z(:);

    overlap=3;

    prepend_x = x(end-overlap:end-1);

    prepend_y = y(end-overlap:end-1);

    prepend_z = z(end-overlap:end-1);

    pospend_x = x(2:overlap+1);

    pospend_y = y(2:overlap+1);

    pospend_z = z(2:overlap+1);

    x = [prepend_x; x; pospend_x];

    y = [prepend_y; y; pospend_y];

    z = [prepend_z; z; pospend_z];

    % SPEED OF CURVE

    dx = gradient(x);

    dy = gradient(y);

    dz = gradient(z);

    dr = [dx dy dz];

    ddx = gradient(dx);

    ddy = gradient(dy);

    ddz = gradient(dz);

    ddr = [ddx ddy ddz];

    % TANGENT

    T = dr./mag(dr,3);

    % DERIVIATIVE OF TANGENT

    dTx = gradient(T(:,1));

    dTy = gradient(T(:,2));

    dTz = gradient(T(:,3));

    dT = [dTx dTy dTz];

    % NORMAL

    N = dT./mag(dT,3);

    % BINORMAL

    B = cross(T,N);

    dBx = gradient(B(:,1));

    dBy = gradient(B(:,2));

    dBz = gradient(B(:,3));

    dB = [dBx dBy dBz];

    % CURVATURE

    % k = mag(dT,1);

    k = mag(cross(dr,ddr),1)./((mag(dr,1)).^3);

    % TORSION

    tor = -dot(dB,N,2);

 

    T = T(overlap+1:end-overlap, :);

    N = N(overlap+1:end-overlap, :);

    B = B(overlap+1:end-overlap, :);

    k = k(overlap+1:end-overlap);

    tor = tor(overlap+1:end-overlap);

end

 

 

function N = mag(T,n)

    % MAGNATUDE OF A VECTOR (Nx3)

    % M = mag(U)

    N = sum(abs(T).^2,2).^(1/2);

    d = find(N==0);

    N(d) = eps*ones(size(d));

    N = N(:,ones(n,1));

end