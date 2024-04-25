clear all;
close all;
flag = "Torus";

aFlag = true;

N = 100;

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

case "MyKnot"

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

end

plot3(x,y,z)

%% TRY ODE45


% Initial curve and time span
% Curv = initial_curve; % Replace initial_curve with your initial curve
y0 = Curv(:); % Reshape the curve to a column vector
tspan = [0,10]; % Time span for simulation

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
lim=3;
grid on;
xlim([-lim, lim]);
ylim([-lim, lim]);
zlim([-10, lim]);

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

% FOR DEBUGGING SHIT

% % Compute the norm of the vectors
% norms = sqrt(sum(Deriv.^2, 1)); % Calculate norm for each vector
% 
% % Plot the norms
% figure
% plot(norms, 'LineWidth', 2);
% xlabel('Index');
% ylabel('Norm');
% title('Norm of Deriv. Vectors along Curve');
% grid on;
% 
% norms2 = sqrt(sum(edgeVec.^2, 1));
% % Plot the norms of edges
% figure
% plot(norms2, 'LineWidth', 2);
% xlabel('Index');
% ylabel('Norm');
% title('Norm of edges along Curve');
% grid on;

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


