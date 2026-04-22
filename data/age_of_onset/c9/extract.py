import pandas as pd

# Read the Excel file (first sheet by default)
df = pd.read_excel('Murphy_2017_c9.xlsx')

# Save to CSV
df.to_csv('Murphy_2017_c9.csv', index=False)

# Read the CSV file
df = pd.read_csv('Murphy_2017_c9.csv')

# Clean up column names (strip whitespace)
df.columns = df.columns.str.strip()

# Filter for ALS and FTD diagnoses
als_mask = df['Diagnosis'].str.contains('ALS', na=False)
ftd_mask = df['Diagnosis'].str.contains('FTD', na=False)

als_df = df[als_mask]
ftd_df = df[ftd_mask]

# For ALS
als_age_counts = als_df['Age of onset'].value_counts().sort_index()
als_age_freq = als_age_counts / als_age_counts.sum() * 100

# For FTD
ftd_age_counts = ftd_df['Age of onset'].value_counts().sort_index()
ftd_age_freq = ftd_age_counts / ftd_age_counts.sum() * 100

print("ALS onset frequency per age (%):")
print(als_age_freq)

print("\nFTD onset frequency per age (%):")
print(ftd_age_freq)

# Define bins (you can adjust the range as needed)
bins = range(int(df['Age of onset'].min()), int(df['Age of onset'].max()) + 6, 5)
labels = [f"{b}-{b+4}" for b in bins[:-1]]

# ALS binned frequencies
als_df['age_bin'] = pd.cut(als_df['Age of onset'], bins=bins, labels=labels, right=False)
als_bin_counts = als_df['age_bin'].value_counts().sort_index()
als_bin_freq = als_bin_counts / als_bin_counts.sum() * 100

print("\nALS onset frequency per 5-year bin (%):")
print(als_bin_freq)
als_bin_freq_df = als_bin_freq.reset_index()
als_bin_freq_df.columns = ['age_bin', 'frequency_percent']
als_bin_freq_df.to_csv('c9_ALS_onset_bin.csv', index=False)

# FTD binned frequencies
ftd_df['age_bin'] = pd.cut(ftd_df['Age of onset'], bins=bins, labels=labels, right=False)
ftd_bin_counts = ftd_df['age_bin'].value_counts().sort_index()
ftd_bin_freq = ftd_bin_counts / ftd_bin_counts.sum() * 100

print("\nFTD onset frequency per 5-year bin (%):")
print(ftd_bin_freq)
ftd_bin_freq_df = ftd_bin_freq.reset_index()
ftd_bin_freq_df.columns = ['age_bin', 'frequency_percent']
ftd_bin_freq_df.to_csv('c9_FTD_onset_bin.csv', index=False)