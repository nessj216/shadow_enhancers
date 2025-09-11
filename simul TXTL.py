import numpy as np
import matplotlib.pyplot as plt

# Constants
gene_lengths = {'sigma': 200, 'anti-sigma': 450}  # in base pairs
transcription_rate = 20  # in bp/sec
translation_rate = 2  # in proteins per bp/sec

# Simulation parameters
total_time = 1000  # in seconds
time_step = 1  # in seconds

# Initialize arrays to store protein levels
sigma_proteins = np.zeros(total_time)
anti_sigma_proteins = np.zeros(total_time)

# Simulation loop
for t in range(total_time):
    # Calculate the number of transcripts produced up to this time step
    sigma_transcripts = transcription_rate * (t + 1)/gene_lengths['sigma']
    anti_sigma_transcripts = transcription_rate * (t + 1)/gene_lengths['anti-sigma']

    # Calculate the number of proteins produced from transcripts
    sigma_proteins[t] = (sigma_transcripts  * translation_rate * (t + 1)) / gene_lengths['sigma']
    anti_sigma_proteins[t] = (anti_sigma_transcripts  * translation_rate * (t + 1)) / gene_lengths['anti-sigma']

# Plot the results
time_points = np.arange(0, total_time, time_step)
plt.plot(time_points, sigma_proteins, label='Sigma Protein')
plt.plot(time_points, anti_sigma_proteins, label='Anti-Sigma Protein')
plt.xlabel('Time (seconds)')
plt.ylabel('Accumulated Protein Levels')
plt.title('Protein Accumulation Simulation')
plt.legend()
plt.show()
