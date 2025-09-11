import numpy as np
import trimesh

# Load the STL file
input_stl_file = "/Users/jillianness/Downloads/Jillian_Slide.stl"

mesh = trimesh.load(input_stl_file)


# Function to scale the length of the square cutout along the x-axis
def scale_square_cutout_length(vertices, xmin, xmax, ymin, ymax, zmin, zmax, scale_factor=0.5):
    # Identify the cutout vertices within the specified x, y, and z range
    cutout_indices = np.where((vertices[:, 0] > xmin) & (vertices[:, 0] < xmax) &
                              (vertices[:, 1] > ymin) & (vertices[:, 1] < ymax) &
                              (vertices[:, 2] > zmin) & (vertices[:, 2] < zmax))

    # Scale the x-coordinates of the identified vertices
    vertices[cutout_indices, 1] *= scale_factor
    return vertices


# Example cutout range (these would be determined from analyzing the mesh)
xmin, xmax = 10, 20  # Update these values based on your specific cutout
ymin, ymax = 10, 20  # Update these values based on your specific cutout
zmin, zmax = 0, 30  # Update these values based on your specific cutout

# Adjust the x-coordinates of the cutout vertices
mesh.vertices = scale_square_cutout_length(mesh.vertices, xmin, xmax, ymin, ymax, zmin, zmax)



# Save the modified STL file
output_stl_file = "/Users/jillianness/Downloads/Modified_Jillian_Slide.stl"
mesh.export(output_stl_file, file_type='stl')

print(f"Modified STL file saved as {output_stl_file}")
