import pandas as pd
import matplotlib.pyplot as plt

# Load your smoothed cumulative frequency data
als_df = pd.read_csv('c9_ALS_onset_smoothed.csv')
als_df['cumulative_frequency'] = als_df['frequency_percent'].cumsum()

# Load the external cumulative incidence data
risk_df = pd.read_csv('/hpc/hers_en/pbeele/simPed/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt', sep='\t')

# Normalize cum_incidence_c9_ALS to its maximum (i.e., set max to 100%)
max_cum_inc = risk_df['cum_incidence_c9_ALS'].max()
risk_df['cum_incidence_c9_ALS_norm'] = risk_df['cum_incidence_c9_ALS'] / max_cum_inc * 100

# Find the age at which the cumulative incidence reaches its maximum
max_idx = risk_df['cum_incidence_c9_ALS'].idxmax()
max_age = risk_df.loc[max_idx, 'age']
max_val = risk_df.loc[max_idx, 'cum_incidence_c9_ALS_norm']

# Plot both curves
plt.figure(figsize=(10, 6))
plt.plot(als_df['age'], als_df['cumulative_frequency'], label='Smoothed ALS Cumulative Frequency (%)')
plt.plot(risk_df['age'], risk_df['cum_incidence_c9_ALS_norm'], label='Cumulative Incidence c9 ALS (Normalized)', linestyle='--')

# Annotate the maximum point
plt.plot(max_age, max_val, 'ro')
plt.axhline(y=max_val, color='red', linestyle=':', linewidth=1)
plt.axvline(x=max_age, color='red', linestyle=':', linewidth=1)
plt.text(max_age, max_val, f'  Max ({max_val:.1f}% at age {max_age})', color='red', va='bottom')

plt.xlabel('Age')
plt.ylabel('Cumulative Percentage (Normalized)')
plt.title('Cumulative ALS Onset vs. Normalized Cumulative Incidence')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('c9_ALS_cumulative_comparison')
plt.show()
