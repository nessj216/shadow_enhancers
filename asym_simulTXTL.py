import numpy as np
import matplotlib.pyplot as plt

# Parameters
transcription_rate = 10  # bp/sec
translation_rate = 5     # bp/sec

sigma_length = 200      # bp
anti_sigma_length = 450  # bp

time = np.linspace(0, 4 * 3600, 100)  # Time in seconds (up to 4 hours)

# Calculate RNA and protein levels over time
sigma_rna = transcription_rate * time / sigma_length
anti_sigma_rna = transcription_rate * time / anti_sigma_length

sigma_protein = translation_rate * time / sigma_length
anti_sigma_protein = translation_rate * time / anti_sigma_length

# Plot
plt.figure(figsize=(10, 6))

plt.plot(time / 3600, sigma_rna, label='Sigma RNA')
plt.plot(time / 3600, anti_sigma_rna, label='Anti-Sigma RNA')

plt.plot(time / 3600, sigma_protein, label='Sigma Protein')
plt.plot(time / 3600, anti_sigma_protein, label='Anti-Sigma Protein')

plt.xlabel('Time (hours)')
plt.ylabel('Relative Levels')
plt.title('Transcription and Translation over Time')
plt.legend()

plt.show()
