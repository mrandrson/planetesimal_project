import pandas as pd
import numpy as np
import os
import glob

directory = "./"  # Set to your directory path if different

for file in glob.glob(os.path.join(directory, "*.csv")):
    df = pd.read_csv(file)

    # Convert 'particle_id' column by extracting integers
    df['particle_id'] = df['particle_id'].astype(str).str.extract(r'(\d+)').astype(float).astype(int)

    # Save the updated dataframe
    df.to_csv(file, index=False)

print("All files updated successfully.")

