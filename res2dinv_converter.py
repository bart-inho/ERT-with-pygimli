import pandas as pd

def read_data_gimli(df):
    # Define column names
    column_names = ['a', 'b', 'm', 'n', 'err', 'i', 'ip', 'iperr', 'k', 'r', 'rhoa', 'u', 'valid']

    # Read the file starting from line 76
    data = pd.read_csv(df, delimiter='\t', skiprows=76, skipfooter=1, names=column_names)

    # Print the data
    print(data)

    return data


def read_geom_gimli(df):
    # Define column names
    column_names = ['x', 'y', 'z']

    # Read the file starting from line 76
    data = pd.read_csv(df, delimiter='\t', skiprows=2, nrows=72, names=column_names)

    print(data)


def convert_2_res2dinv(df, filename='data_res2dinv/ert_WIL_01.dat'):

    # Create a new dataframe
    new_df = pd.DataFrame(columns=['id', 'a', 'type_a', 'b', 'type_b', 'm', 'type_m', 'n', 'type_n', 'rhoa'])

    new_df['a'] = (df['a']-1) * 2.5
    new_df['b'] = (df['b']-1) * 2.5
    new_df['m'] = (df['m']-1) * 2.5
    new_df['n'] = (df['n']-1) * 2.5

    new_df['rhoa'] = df['rhoa']

    new_df['type_a'] = 0.0
    new_df['type_b'] = 0.0
    new_df['type_m'] = 0.0
    new_df['type_n'] = 0.0

    new_df['id'] = int(4)

    header = """P01_O-E_clean-topo.bin 
         2.50
11
7
Type of measurement (0=app. resistivity,1=resistance)
0
1316
 2
0
"""

    footer = """0
0
0
0
"""

    # Open the file in write mode
    with open(filename, 'w') as f:
        # Write the header
        f.write(header)
        
        # Write each row to the file with the desired spacing
        for _, row in new_df.iterrows():
            line = f"{int(row['id']):1}  {row['a']:2.2f}  {row['type_a']:.2f} {row['b']:2.2f}  {row['type_b']:.2f} {row['m']:2.2f}  {row['type_m']:.2f} {row['n']:2.2f}  {row['type_n']:.2f} {row['rhoa']:2.5f}\n"
            f.write(line)

        
        f.write(footer)


    print(new_df)
    
    return new_df

# def read_res2Dinv_processed(path):
#     import matplotlib.pyplot as plt

#     # Read the data, skipping the first 5 rows and using whitespace as the delimiter
#     df = pd.read_csv(path, skiprows=2473, delim_whitespace=True, nrows=1136, names=['X', 'Depth', 'Resistivity', 'Relative sensitivity', 'Smoothed sensitivity', 'Uncertainty'])

#     # Pivot the DataFrame to create a grid of resistivity values
#     pivot_df = df.pivot(index='Depth', columns='X', values='Resistivity')

#     # Plot the colormap using the 'jet' cmap
#     plt.figure(figsize=(10, 8))
#     plt.imshow(pivot_df, cmap='jet', origin='lower')
#     plt.colorbar(label='Resistivity')
#     plt.title('Resistivity Colormap')
#     plt.xlabel('X')
#     plt.ylabel('Depth')
#     plt.show()