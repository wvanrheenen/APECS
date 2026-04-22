import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Read the input file (tab-delimited, adjust sep=',' if comma-delimited)
df = pd.read_csv('Nona_polygenicALS.txt', sep='\t')  # or sep=',' if needed

# Sort by age to ensure correct interpolation
df = df.sort_values('age')

# Interpolation
interp_func = interp1d(df['age'], df['frequency'], kind='cubic', bounds_error=False, fill_value=0)

# New ages (every integer from 12 to 92)
ages_new = np.arange(15, 93)
freqs_new = interp_func(ages_new)

# Set negative interpolated frequencies to zero
freqs_new = np.clip(freqs_new, 0, None)

# Normalize so total frequency sums to 100%
freqs_new = freqs_new / freqs_new.sum() * 100

# Save to CSV
df_interp = pd.DataFrame({'age': ages_new, 'frequency': freqs_new})
df_interp.to_csv('als_onset_polygenic.csv', index=False)

df = pd.read_csv('als_onset_polygenic.csv')
total_frequency = df['frequency'].sum()
print(total_frequency)

# Calculate the cumulative frequency
df['cumulative_frequency'] = df['frequency'].cumsum()

# Plot cumulative frequency
plt.figure(figsize=(8, 5))
plt.plot(df['age'], df['cumulative_frequency'], label='Polygenic ALS (Cumulative)')
plt.xlabel('Age')
plt.ylabel('Cumulative Frequency (%)')
plt.title('Cumulative Frequency of ALS Onset (Polygenic)')
plt.grid(True)
plt.tight_layout()
plt.savefig('polygenic_ALS_cumulative.png')
plt.show()

# Load interpolated data
interp_df = pd.read_csv('als_onset_polygenic.csv')

plt.figure(figsize=(10,5))
plt.bar(interp_df['age'], interp_df['frequency'], width=1, alpha=0.6, label='Interpolated (per age)')
plt.xlabel('Age')
plt.ylabel('Frequency (%)')
plt.title('ALS Onset Frequency Distribution (Age-Specific)')
plt.legend()
plt.tight_layout()
plt.savefig('polygenic_histo.png')
plt.show()