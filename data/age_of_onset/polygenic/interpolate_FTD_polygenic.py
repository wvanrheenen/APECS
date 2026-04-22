import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Read the FTD input file (comma-delimited as provided)
df = pd.read_csv('Polygenic_FTD_incidence.csv')  # Assumes file with columns: Age, corrected_incidence

print("Actual column names:", df.columns.tolist())
print("\nFirst 3 rows:\n", df.head(3))
print("\nData types:\n", df.dtypes)
print("\nShape:", df.shape)

# Sort by age to ensure correct interpolation
df = df.sort_values('age')

# Interpolation (cubic spline for smooth curve)
interp_func = interp1d(df['age'], df['frequency'], kind='cubic', bounds_error=False, fill_value=0)

# New ages (every integer from 35 to 100 for FTD range)
ages_new = np.arange(35, 101)
freqs_new = interp_func(ages_new)

# Set negative interpolated frequencies to zero
freqs_new = np.clip(freqs_new, 0, None)

# Normalize so total frequency sums to 100%
freqs_new = freqs_new / freqs_new.sum() * 100

# Save interpolated frequencies to CSV
df_interp = pd.DataFrame({'age': ages_new, 'frequency': freqs_new})
df_interp.to_csv('ftd_onset_polygenic.csv', index=False)

# Verify total sums to ~100%
df = pd.read_csv('ftd_onset_polygenic.csv')
total_frequency = df['frequency'].sum()
print(f"Total frequency: {total_frequency:.2f}%")

# Calculate cumulative frequency (now properly represents smoothed polygenic FTD cumulative incidence)
df['cumulative_frequency'] = df['frequency'].cumsum()

# Plot cumulative incidence (like Nona Fig 1)
plt.figure(figsize=(8, 5))
plt.plot(df['age'], df['cumulative_frequency'], 'b-', linewidth=2, label='Polygenic FTD (Cumulative Incidence)')
plt.xlabel('Age at Onset (years)')
plt.ylabel('Cumulative Incidence (%)')
plt.title('Cumulative Incidence of FTD Onset (Age-Corrected, Polygenic)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig('polygenic_FTD_cumulative.png', dpi=300, bbox_inches='tight')
plt.show()

# Plot age-specific frequency distribution (histogram)
plt.figure(figsize=(10, 5))
plt.bar(df['age'], df['frequency'], width=1, alpha=0.7, color='lightblue', edgecolor='navy', linewidth=0.5)
plt.xlabel('Age at Onset (years)')
plt.ylabel('Frequency (%)')
plt.title('FTD Onset Frequency Distribution (Interpolated, Age-Corrected)')
plt.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig('polygenic_FTD_histo.png', dpi=300, bbox_inches='tight')
plt.show()
