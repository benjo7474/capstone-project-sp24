import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from findiff import FinDiff

# Define the circle and motion parameters
R = 1  # Radius of the circle
num_points = 100  # Number of points along the circle
dt = 0.05  # Time step
total_time = 2.0  # Total time for the simulation
time_steps = int(total_time / dt)
curvature = 1 / R  # Constant curvature

# Parametrize the circle
t = np.linspace(0, 2 * np.pi, num_points)
x = R * np.cos(t)
y = R * np.sin(t)
z = np.zeros(num_points)

# Define differential operator for derivative calculation
dx = FinDiff(0, t[1] - t[0], 1)

# Compute the Frenet frame
T = np.vstack((dx(x), dx(y), np.zeros(num_points)))  # Using findiff for derivatives
norm_T = np.linalg.norm(T, axis=0)
T = T / norm_T

# Calculate dT/ds using np.gradient since it involves derivatives of unit vectors
dT = np.gradient(T, axis=1)
N = np.cross(np.cross(T, dT, axis=0), T, axis=0)
norm_N = np.linalg.norm(N, axis=0)
N = N / norm_N
B = np.cross(T, N, axis=0)
norm_B = np.linalg.norm(B, axis=0)
B = B / norm_B

# Initialize positions
X = np.vstack((x, y, z))
positions = [X]

def update_position(X, B, curvature, dt):
    return X + curvature * B * dt

for _ in range(time_steps):
    X = update_position(X, B, curvature, dt)
    positions.append(X)

# Animation setup
fig, ax = plt.subplots(subplot_kw={'projection': '3d'}, figsize=(10, 8))
line, = ax.plot([], [], [], 'b-')

def init():
    line.set_data([], [])
    line.set_3d_properties([])
    return line,

def update(frame):
    pos = positions[frame]
    line.set_data(pos[0, :], pos[1, :])
    line.set_3d_properties(pos[2, :])
    return line,

ani = animation.FuncAnimation(fig, update, frames=len(positions),
                              init_func=init, blit=True, interval=50)

ax.set_xlim([-1.5, 1.5])
ax.set_ylim([-1.5, 1.5])
ax.set_zlim([-1.5, 1.5])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Time Evolution of a Vortex Filament (Circle)')

# Save or display the animation
ani.save('vortex_filament_evolution.mp4', writer='ffmpeg', dpi=200)
# plt.show() to view inline if running interactively

plt.close()  # Prevents duplicate display in notebooks
