import pandas as pd

# If the file has no header, specify header=None and assign column names
df_check = pd.read_csv('c9_FTD_check_Moore.csv', header=None, names=['age', 'cum_freq'])
df_check.to_csv('c9_FTD_check_Moore.csv', index=False)  # Overwrite with header

import matplotlib.pyplot as plt

# Load the smoothed data
df_smooth = pd.read_csv('c9_FTD_onset_smoothed.csv')
df_smooth['cumulative_frequency'] = df_smooth['frequency_percent'].cumsum()

# Load the check Moore data (now with headers)
df_check = pd.read_csv('c9_FTD_check_Moore.csv')

# Plot both cumulative frequencies
plt.figure(figsize=(10, 6))
plt.plot(df_smooth['age'], df_smooth['cumulative_frequency'], label='Smoothed C9 FTD (cumulative)', linewidth=2)
plt.plot(df_check['age'], df_check['cum_freq'], label='Check Moore C9 FTD (cumulative)', linewidth=2, linestyle='--')
plt.xlabel('Age')
plt.ylabel('Cumulative Frequency (%)')
plt.title('Cumulative Frequency: Smoothed vs. Check Moore C9 FTD')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('c9_FTD_cumulative_comparison.png')
plt.show()