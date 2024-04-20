import numpy as np
import matplotlib.pyplot as plt

# Define the parameters for the circle
R = 1  # Radius of the circle
num_points = 100  # Number of points along the circle

# Parametrize the circle
t = np.linspace(0, 2 * np.pi, num_points)  # Parameter t from 0 to 2*pi
x = R * np.cos(t)  # x-coordinate
y = R * np.sin(t)  # y-coordinate
z = np.zeros(num_points)  # z-coordinate (0 since the circle is in the XY-plane)

# Compute the derivatives for the tangent vector
dx = -R * np.sin(t)
dy = R * np.cos(t)
dz = np.zeros(num_points)

# Normalize the tangent vector (T)
T = np.vstack((dx, dy, dz))
norm_T = np.linalg.norm(T, axis=0)
T = T / norm_T  # Normalized tangent vector

# Compute the derivative of T to find the normal vector (N)
dT = np.gradient(T, axis=1)  # Numerical derivative of T
N = np.cross(np.cross(T, dT, axis=0), T, axis=0)  # Double cross product to ensure orthogonality and direction
norm_N = np.linalg.norm(N, axis=0)
N = N / norm_N  # Normalized normal vector

# Compute the binormal vector (B) as the cross product of T and N
B = np.cross(T, N, axis=0)
norm_B = np.linalg.norm(B, axis=0)
B = B / norm_B  # Normalized binormal vector

# Visualization of the Frenet Frame along the circle
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, label='Circle', color='blue')  # Plot the circle

# Plot Frenet frames at selected points
indices = [10, 30, 50, 70, 90]  # Indices of points to plot frames
scale = 0.2  # Scaling factor for the frame vectors
for i in indices:
    ax.quiver(x[i], y[i], z[i], T[0, i], T[1, i], T[2, i], color='red', length=scale, label='Tangent' if i == indices[0] else "")
    ax.quiver(x[i], y[i], z[i], N[0, i], N[1, i], N[2, i], color='green', length=scale, label='Normal' if i == indices[0] else "")
    ax.quiver(x[i], y[i], z[i], B[0, i], B[1, i], B[2, i], color='black', length=scale, label='Binormal' if i == indices[0] else "")

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.title('Circle with Frenet Frame Vectors')
plt.show()
