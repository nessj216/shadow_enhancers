import numpy as np
import matplotlib.pyplot as plt

# Constants and parameters
x_positions = np.linspace(0, 50, 100)  # X positions along the embryo length (0 to 10 units)
fixed_y_position = 1  # Fixed Y position
correct_ap_angle = 0  # Correct AP angle in radians
error_angles = [0.1, 0.2, 0.3]  # Different levels of angle error in radians

# Function to calculate AP position error
def calculate_error(x_positions, y_position, angle_error):
    distances = np.sqrt(x_positions**2 + y_position**2)
    angles = np.arctan2(y_position, x_positions)
    return distances * np.abs(np.cos(angles + correct_ap_angle) - np.cos(angles + correct_ap_angle + angle_error))

# Plotting
plt.figure(figsize=(10, 6))

for error_angle in error_angles:
    errors = calculate_error(x_positions, fixed_y_position, error_angle)
    plt.plot(x_positions, errors, label=f'Angle error = {error_angle} radians')

plt.xlabel('X Position along Embryo Length')
plt.ylabel('Estimated Error in AP Position')
plt.title('Error in AP Position as a Function of X Position')
plt.legend()
plt.grid(True)
plt.show()
