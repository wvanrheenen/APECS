import pandas as pd

# Load TSV file
df = pd.read_csv('lifetime_dementia_risk.txt', sep='\t')

# Calculate annual frequency (difference in cumulative incidence)
df['frequency'] = df['cumulative_incidence_perc'].diff().fillna(df['cumulative_incidence_perc'])

# Normalize frequencies so they sum to 100%
df['frequency'] = df['frequency'] / df['frequency'].sum() * 100

# Save as CSV with age and normalized frequency columns
df[['age', 'frequency']].to_csv('polygenic_dementia.csv', index=False)