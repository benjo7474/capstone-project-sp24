import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from findiff import FinDiff

# Define the helix and motion parameters
radius = 1  # Radius of the main helix
small_radius = 0.2  # Radius of the high-frequency corkscrew
pitch = 2  # Vertical distance for a single complete loop of the helix
high_freq = 10  # Frequency multiplier for the high-frequency corkscrew
num_points = 1000  # Number of points along the helix
dt = 0.05  # Time step
total_time = 10.0  # Total time for the simulation
time_steps = int(total_time / dt)

# Parametrize the helix
t = np.linspace(0, 4 * np.pi, num_points)  # 2 loops
x = radius * np.cos(t) + small_radius * np.cos(high_freq * t)
y = radius * np.sin(t) + small_radius * np.sin(high_freq * t)
z = pitch * t / (2 * np.pi)

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
ax.set_title('The Cash Special: Helix Corkscrew Evolution')

# Save or display the animation
ani.save('the_cash_special.mp4', writer='ffmpeg', dpi=200)
# plt.show() to view inline if running interactively

plt.close()  # Prevents duplicate display in notebooks
