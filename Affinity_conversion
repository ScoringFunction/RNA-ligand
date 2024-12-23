import pandas as pd
import numpy as np

def convert_affinity_to_kcal(affinity, unit):
    """
    Convert binding affinity from nM, uM, mM, or M to kcal/mol.
    
    Parameters:
    affinity (float): The binding affinity value.
    unit (str): The unit of the binding affinity ('nM', 'uM', 'mM', 'M').
    
    Returns:
    float: The binding affinity in kcal/mol (negative value).
    """
    # Constants
    R = 1.987  # cal/(mol*K)
    T = 298.15  # K (25°C)
    
    # Convert affinity to Molarity (M)
    if unit == 'nM':
        affinity_molar = affinity * 1e-9
    elif unit == 'uM':
        affinity_molar = affinity * 1e-6
    elif unit == 'mM':
        affinity_molar = affinity * 1e-3
    elif unit == 'M':
        affinity_molar = affinity
    else:
        raise ValueError("Invalid unit. Please use 'nM', 'uM', 'mM', or 'M'.")
    
    # Calculate binding affinity in kcal/mol
    delta_G = R * T * np.log(affinity_molar) / 1000  # Convert cal to kcal
    
    return -abs(delta_G)

# Load the binding affinities from the provided Excel file
binding_affinities_df = pd.read_excel('binding_affinities.xlsx')

# Extract the PDB File and binding affinities
pdb_files = binding_affinities_df['PDB File']
binding_affinities = binding_affinities_df['binding_affinity']

# Initialize a list to store the converted affinities
converted_affinities = []

# Convert each binding affinity to kcal/mol
for affinity in binding_affinities:
    # Ensure the affinity is a string before processing
    affinity_str = str(affinity).strip()
    
    # Extract the numeric value and unit from the affinity string
    if 'nM' in affinity_str:
        value = float(affinity_str.replace('nM', '').strip())
        unit = 'nM'
    elif 'uM' in affinity_str:
        value = float(affinity_str.replace('uM', '').strip())
        unit = 'uM'
    elif 'mM' in affinity_str:
        value = float(affinity_str.replace('mM', '').strip())
        unit = 'mM'
    elif 'M' in affinity_str:
        value = float(affinity_str.replace('M', '').strip())
        unit = 'M'
    else:
        # Handle scientific notation and other edge cases
        try:
            value = float(affinity_str)
            unit = 'M'
        except ValueError:
            raise ValueError(f"Invalid unit in binding affinity: {affinity_str}. Please use 'nM', 'uM', 'mM', or 'M'.")
    
    # Convert the affinity to kcal/mol
    converted_affinity = convert_affinity_to_kcal(value, unit)
    
    # Append the converted affinity to the list
    converted_affinities.append(converted_affinity)

# Create a new DataFrame with the PDB File and converted affinities
converted_affinities_df = pd.DataFrame({
    'PDB File': pdb_files,
    'Binding Affinity (kcal/mol)': converted_affinities
})

# Save the converted affinities to a new Excel file
converted_affinities_df.to_excel('Converted_Binding_Affinities.xlsx', index=False)

print("Converted binding affinities saved to Converted_Binding_Affinities.xlsx")
