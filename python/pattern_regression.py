



pattern_csv = '../pvcsummary_bl_patterns_and_change.csv'

df = pd.read_csv(pattern_csv)
df.set_index('rid',inplace=True)