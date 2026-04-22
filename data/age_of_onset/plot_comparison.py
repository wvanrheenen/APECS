import pandas as pd
import matplotlib.pyplot as plt

# Load data 
als_poly_df = pd.read_csv('polygenic/als_onset_polygenic.csv')
ftd_poly_df = pd.read_csv('polygenic/ftd_onset_polygenic.csv')
dem_poly_df = pd.read_csv('polygenic/polygenic_dementia.csv')  # ADDED
c9_als_df = pd.read_csv('c9/c9_ALS_onset_smoothed.csv')
c9_ftd_df = pd.read_csv('c9/c9_FTD_onset_smoothed.csv')
c9_dementia_df = pd.read_csv('c9/c9_dementia_smoothed.csv')

# Calculate cumulative frequencies
dfs = [als_poly_df, ftd_poly_df, dem_poly_df, c9_als_df, c9_ftd_df, c9_dementia_df]  # UPDATED
cum_cols = ['frequency', 'frequency', 'frequency', 'frequency_percent', 'frequency_percent', 'frequency']
for df, col in zip(dfs, cum_cols):
    df['cumulative_frequency'] = df[col].cumsum()

# Plot all cumulative frequencies (ADDED dem_poly_df)
plt.figure(figsize=(12, 7))
plt.plot(als_poly_df['age'], als_poly_df['cumulative_frequency'], label='Polygenic ALS', linewidth=2, color='skyblue', linestyle='-')
plt.plot(ftd_poly_df['age'], ftd_poly_df['cumulative_frequency'], label='Polygenic FTD', linewidth=2)
plt.plot(dem_poly_df['age'], dem_poly_df['cumulative_frequency'], label='Polygenic Dementia',  # ADDED
         linewidth=2, linestyle='-', color='lightgreen')
plt.plot(c9_als_df['age'], c9_als_df['cumulative_frequency'], label='C9 ALS', linewidth=2, color='darkblue', linestyle='--')
plt.plot(c9_ftd_df['age'], c9_ftd_df['cumulative_frequency'], label='C9 FTD', linewidth=2, color='red', linestyle='--')
plt.plot(c9_dementia_df['age'], c9_dementia_df['cumulative_frequency'], 
         label='C9 Dementia', linewidth=2, linestyle='--', color='darkgreen')

plt.xlabel('Age', fontsize=12)
plt.ylabel('Cumulative Frequency (%)', fontsize=12)
plt.title('Neurodegenerative Disease Onset Comparison', fontsize=14)
plt.legend(fontsize=10)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('cumulative_comparison.png', dpi=300)
plt.show()

# Prepare merged CSV (ADDED dem_poly_df)
merged = (
    als_poly_df[['age', 'frequency']].rename(columns={'frequency': 'Polygenic_ALS'})
    .merge(ftd_poly_df[['age', 'frequency']].rename(columns={'frequency': 'Polygenic_FTD'}), 
            on='age', how='outer')
    .merge(dem_poly_df[['age', 'frequency']].rename(columns={'frequency': 'Polygenic_Dementia'}),  # ADDED
            on='age', how='outer')
    .merge(c9_als_df[['age', 'frequency_percent']].rename(columns={'frequency_percent': 'C9_ALS'}), 
            on='age', how='outer')
    .merge(c9_ftd_df[['age', 'frequency_percent']].rename(columns={'frequency_percent': 'C9_FTD'}), 
            on='age', how='outer')
    .merge(c9_dementia_df[['age', 'frequency']].rename(columns={'frequency': 'C9_Dementia'}), 
            on='age', how='outer')
)

# Clean and save
merged.sort_values('age').fillna(0).to_csv('age_of_onset.csv', index=False)
