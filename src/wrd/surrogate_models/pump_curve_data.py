import pandas as pd
import json

# TODO: Convert to functions

# # For RO pumps (no efficiency curve)
# # Load the JSON file as a dictionary first
# with open(r'C:\Users\rchurchi\flex_desal\src\wrd\surrogate_models\RO_IS_pump_curves.json', 'r') as f:
#     data = json.load(f)

# datasets = data['datasetColl']

# # Convert specific datasets to DataFrames
# eff_data = pd.DataFrame()
# head_data = pd.DataFrame()
# for dataset in datasets:
#     name = dataset['name']
#     df = pd.DataFrame(dataset['data'])
#     try:
#         eff = int(name[:2])  # Extract first 2 characters and convert to int
#         df['efficiency'] = eff  # add efficiency column with every value the same
#         df['Flow (gpm)'] = df['value'].apply(lambda x: x[0])
#         df['Head (ft)'] = df['value'].apply(lambda x: x[1])
#         df.drop(columns=['x','y','value'], inplace=True)
#         eff_data = pd.concat([eff_data, df], ignore_index=True)
#     except:
#         df['Flow (gpm)'] = df['value'].apply(lambda x: x[0])
#         df[name+' (ft)'] = df['value'].apply(lambda x: x[1])
#         df.drop(columns=['x','y','value'], inplace=True)
#         head_data = pd.concat([head_data, df], ignore_index=True)


# print(eff_data.head())
# eff_data.to_csv(r'C:\Users\rchurchi\flex_desal\src\wrd\surrogate_models\RO_IS_pump_eff_curve_data.csv', index=False)

# print(head_data.head())
# head_data.to_csv(r'C:\Users\rchurchi\flex_desal\src\wrd\surrogate_models\RO_IS_pump_head_curves_data.csv', index=False)

# For UF pump (with efficiency curve)
with open(r'C:\Users\rchurchi\flex_desal\src\wrd\surrogate_models\uf_pump_curves.json', 'r') as f:
    data = json.load(f)

datasets = data['datasetColl']
# Convert specific datasets to DataFrames
eff_data = pd.DataFrame()
head_data = pd.DataFrame()
for dataset in datasets:
    name = dataset['name']
    print(name)
    df = pd.DataFrame(dataset['data'])
    if name == 'Efficiency':
        df['Flow (gpm)'] = df['value'].apply(lambda x: x[0])
        df['Efficiency (%)'] = df['value'].apply(lambda x: x[1])
        df.drop(columns=['x','y','value'], inplace=True)
        eff_data = pd.concat([eff_data, df], ignore_index=True)        
    else:
        df['Flow (gpm)'] = df['value'].apply(lambda x: x[0])
        df[name + ' (ft)'] = df['value'].apply(lambda x: x[1])
        df.drop(columns=['x','y','value'], inplace=True)
        head_data = pd.concat([head_data, df], ignore_index=True)
print(eff_data.head())
eff_data.to_csv(r'C:\Users\rchurchi\flex_desal\src\wrd\surrogate_models\RO_UF_pump_eff_curve_data.csv', index=False)
print(head_data.head())
head_data.to_csv(r'C:\Users\rchurchi\flex_desal\src\wrd\surrogate_models\RO_UF_pump_head_curves_data.csv', index=False)