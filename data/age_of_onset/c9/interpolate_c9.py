import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

def smooth_binned_onset(input_csv, output_csv, age_start=22, age_end=97, plot_title=None):
    # Read the binned data
    df = pd.read_csv(input_csv)
    
    # Calculate bin midpoints
    df['age_mid'] = df['age_bin'].str.extract(r'(\d+)-(\d+)').astype(int).mean(axis=1)
    
    # Normalize frequencies to proportions
    df['frequency_prop'] = df['frequency_percent'] / df['frequency_percent'].sum()
    
    # Interpolation range
    ages_new = np.arange(age_start, age_end + 1)
    
    # Spline interpolation (quadratic for stability)
    spline = make_interp_spline(df['age_mid'], df['frequency_prop'], k=2)
    freqs_new = spline(ages_new)
    freqs_new = np.clip(freqs_new, 0, None)
    freqs_new /= freqs_new.sum()  # Normalize to sum to 1
    freqs_new_percent = freqs_new * 100
    
    # Save to CSV
    pd.DataFrame({'age': ages_new, 'frequency_percent': freqs_new_percent}).to_csv(output_csv, index=False)

# Example usage:
smooth_binned_onset('c9_ALS_onset_bin.csv', 'c9_ALS_onset_smoothed.csv', age_start=22, age_end=87, plot_title='C9 ALS Onset')
smooth_binned_onset('c9_FTD_onset_bin.csv', 'c9_FTD_onset_smoothed.csv', age_start=22, age_end=87, plot_title='C9 FTD Onset')


# ALS
als_df = pd.read_csv('c9_ALS_onset_smoothed.csv')
als_df['cumulative_frequency'] = als_df['frequency_percent'].cumsum()

plt.figure(figsize=(8, 5))
plt.plot(als_df['age'], als_df['cumulative_frequency'], label='C9 ALS (Cumulative)')
plt.xlabel('Age')
plt.ylabel('Cumulative Frequency (%)')
plt.title('Cumulative Frequency of ALS Onset (C9)')
plt.grid(True)
plt.tight_layout()
plt.savefig('c9_ALS_cumulative.png')
plt.show()

# FTD
ftd_df = pd.read_csv('c9_FTD_onset_smoothed.csv')
ftd_df['cumulative_frequency'] = ftd_df['frequency_percent'].cumsum()

plt.figure(figsize=(8, 5))
plt.plot(ftd_df['age'], ftd_df['cumulative_frequency'], label='C9 FTD (Cumulative)', color='orange')
plt.xlabel('Age')
plt.ylabel('Cumulative Frequency (%)')
plt.title('Cumulative Frequency of FTD Onset (C9)')
plt.grid(True)
plt.tight_layout()
plt.savefig('c9_FTD_cumulative.png')
plt.show()