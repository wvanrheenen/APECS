import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

# Step 1: Read the CSV and add headers
df = pd.read_csv('c9_dementia.csv', header=None, names=['age', 'cum_frequency'])

# Step 2: Round ages and group by rounded age, averaging cum_frequency for duplicates
df['age_rounded'] = df['age'].round().astype(int)
df_grouped = df.groupby('age_rounded')['cum_frequency'].mean().reset_index()

# Step 3: Add start and end points
df_grouped = pd.concat([
    pd.DataFrame({'age_rounded': [35], 'cum_frequency': [0.0]}),
    df_grouped,
    pd.DataFrame({'age_rounded': [85], 'cum_frequency': [1.0]})
], ignore_index=True).sort_values('age_rounded')

# Step 4: Ensure monotonicity before interpolation
df_grouped['cum_frequency'] = np.maximum.accumulate(df_grouped['cum_frequency'])

# Step 5: Interpolate with PCHIP for all integer ages between 35 and 85
age_grid = np.arange(35, 86)
pchip = PchipInterpolator(df_grouped['age_rounded'], df_grouped['cum_frequency'])
cum_freq_grid = pchip(age_grid)
cum_freq_grid = np.clip(cum_freq_grid, 0, 1)

# Step 6: Calculate frequency per age (difference between cumulative frequencies)
frequency_grid = np.diff(cum_freq_grid, prepend=0)  # prepend 0 for the first age

# Step 7: Build result DataFrame
result_df = pd.DataFrame({
    'age': age_grid,
    'cum_frequency': cum_freq_grid,
    'frequency': frequency_grid
})
result_df['cum_frequency'] = result_df['cum_frequency'] * 100
result_df['frequency'] = result_df['frequency'] * 100

# Step 8: Plot (optional)
plt.figure(figsize=(8, 5))
plt.plot(result_df['age'], result_df['cum_frequency'], '-', label='Cumulative Frequency')
plt.scatter(df['age'], df['cum_frequency'], color='red', s=15, label='Original Data Points', zorder=5)
plt.xlabel('Age')
plt.ylabel('Cumulative Frequency')
plt.title('Monotonic Cumulative Incidence Curve (C9 Dementia)')
plt.legend()
plt.savefig('c9_dementia_cumulative.png')
plt.show()

# Step 9: Save to CSV
result_df.to_csv('c9_dementia_smoothed.csv', index=False)