import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from findiff import FinDiff

# Define the Torus Knot and motion parameters
num_points = 100000  # Number of points along the Torus Knot
dt = 0.1  # Time step
total_time = 10.0  # Total time for the simulation
time_steps = int(total_time / dt)

# Parametrize the Torus Knot
t = np.linspace(0, 4 * np.pi, num_points)  # Four loops
x = 1/3 * (np.cos(np.pi * 2 * t) + 2 * np.cos(2 * np.pi * 2 * t))
y = 1/3 * (np.sin(np.pi * 2 * t) - 2 * np.sin(2 * np.pi * 2 * t))
z = 1/3 * (2 * np.sin(np.pi * 2 * t))

# Define differential operator for derivative calculation
dx = FinDiff(0, t[1] - t[0], 1)

# Compute the Frenet frame
T = np.vstack((dx(x), dx(y), dx(z)))  # Using findiff for derivatives
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

def update_position(X, B, dt):
    return X + B * dt

for _ in range(time_steps):
    X = update_position(X, B, dt)
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

# Calculate maximum extents
all_positions = np.hstack(positions)
x_limits = [np.min(all_positions[0, :]), np.max(all_positions[0, :])]
y_limits = [np.min(all_positions[1, :]), np.max(all_positions[1, :])]
z_limits = [np.min(all_positions[2, :]), np.max(all_positions[2, :])]

ax.set_xlim(x_limits)
ax.set_ylim(y_limits)
ax.set_zlim(z_limits)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Time Evolution of a Helix')

# Save or display the animation
ani.save('torus_evolution.mp4', writer='ffmpeg', dpi=200)
# plt.show() to view inline if running interactively

plt.close()  # Prevents duplicate display in notebooks
