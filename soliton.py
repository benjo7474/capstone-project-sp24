import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from findiff import FinDiff


def gaussian_kernel(size, sigma=1.0):
    """ Returns a 1D Gaussian kernel """
    size = int(size) // 2
    x = np.linspace(-size, size, 2*size + 1)
    g = np.exp(-(x**2 / (2 * sigma**2)))
    return g / g.sum()

# Parameters
sigma = 100.0  # Standard deviation of the Gaussian kernel
kernel_size = int(6*sigma + 1)  # Size of the kernel

# Generate Gaussian kernel
kernel = gaussian_kernel(kernel_size, sigma)

# Define the helix and motion parameters
radius = 10  # Radius of the main helix
small_radius = 0.2  # Radius of the high-frequency corkscrew
pitch = 15  # Vertical distance for a single complete loop of the helix
high_freq = 3  # Frequency multiplier for the high-frequency corkscrew
num_points = 10000  # Number of points along the helix
dt = 0.05  # Time step
total_time = 10.0  # Total time for the simulation
time_steps = int(total_time / dt)

# Parametrize the helix
t = np.linspace(0 * np.pi, 2.5 * np.pi, num_points)  # 2 loops


gaussian = 10*np.exp(-(t-1.5*np.pi)**2)

x = t+gaussian*(np.cos(high_freq*t)+np.sin(high_freq*t))
y = t+gaussian*(np.cos(high_freq*t)-np.sin(high_freq*t))
z = t+gaussian*(-np.cos(high_freq*t))#small_radius*np.sin(high_freq * t)

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

# Use fewer points for arrows to avoid clutter
arrow_indices = np.linspace(0, num_points - 1, 300, dtype=int)  # Adjust number for more or fewer arrows

for _ in range(time_steps):
    X = update_position(X, B, dt)
    positions.append(X)

# Animation setup
fig, ax = plt.subplots(subplot_kw={'projection': '3d'}, figsize=(10, 8))
line, = ax.plot([], [], [], 'b-')
quiver = None

def init():
    line.set_data([], [])
    line.set_3d_properties([])
    global quiver
    if quiver:
        quiver.remove()
    quiver = ax.quiver([], [], [], [], [], [], color='r', length=0.5, normalize=True)
    return line, quiver

def update(frame):
    pos = positions[frame]
    line.set_data(pos[0, :], pos[1, :])
    line.set_3d_properties(pos[2, :])
    # Update arrows for binormal vectors
    global quiver
    if quiver:
        quiver.remove()
    bx, by, bz = B[0, arrow_indices], B[1, arrow_indices], B[2, arrow_indices]
    quiver = ax.quiver(pos[0, arrow_indices], pos[1, arrow_indices], pos[2, arrow_indices],
                       bx, by, bz, color='r', length=0.5, normalize=True)
    return line, quiver

ani = animation.FuncAnimation(fig, update, frames=len(positions),
                              init_func=init, blit=True, interval=50)

# Calculate maximum extents
all_positions = np.hstack(positions)
x_limits = [np.min(all_positions[0, :]), np.max(all_positions[0, :])]
y_limits = [np.min(all_positions[1, :]), np.max(all_positions[1, :])]
z_limits = [np.min(all_positions[2, :]), np.max(all_positions[2, :])]

print(x_limits)

ax.set_xlim(x_limits)
ax.set_ylim(y_limits)
ax.set_zlim(z_limits)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Soliton with Labeled Binormal')

# Save or display the animation
ani.save('soliton.mp4', writer='ffmpeg', dpi=200)
# plt.show() to view inline if running interactively

plt.close()  # Prevents duplicate display in notebooks
