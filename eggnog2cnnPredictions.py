# Import packages:
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
import argparse

# Add menu for command line usage:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run eggnog to CNN result parameters')
    parser.add_argument('-file', action='store', dest='file', required=True,
                        help='eggnog output file')
    parser.add_argument('-o', action='store', dest='output_file',
                        help='Name of output file ext .txt', required=True)
    parser.add_argument('-m', action='store', dest='model',
                        help='Name of CNN model', required=True)

options = parser.parse_args()

#Read in eggnog file:
eggnog = pd.read_csv(options.file, sep='\\t', comment='#',names=['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name', 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'] )
# Subset table only selecting column 1 and 5
eggnog_cogs = eggnog[['query', 'eggNOG_OGs']]

eggnog_cogs['query'] = eggnog_cogs['query'].str.split('_').str[0]
eggnog_cogs['eggNOG_OGs'] = eggnog_cogs['eggNOG_OGs'].apply(lambda s: s.split('@')[0])

eggnog_cogs = eggnog_cogs.rename(columns={'eggNOG_OGs': 'cogs'})

# Add column with 'your_genome' to avoid weirdness with contig ids within the genomes
eggnog_cogs['genome_id'] = 'your_genome'

grouped = eggnog_cogs.groupby(["genome_id", "cogs"]).size().reset_index(name="cog_freq")

table = pd.pivot_table(grouped, values=['cog_freq'], index=['genome_id'],
                    columns=['cogs'], fill_value=0)

table = table.reset_index()
table.columns = table.columns.droplevel(0)
table.rename(columns = {list(table)[0]: 'genome_id'}, inplace = True)
# Get a unique list of the genomes:
genome_ids = table['genome_id'].values

# Read in the gene families which the CNN was trained on:
with open('genefamilies.txt', 'r') as file:
    desired_columns = [line.strip() for line in file.readlines()]

# Ensure that the column names in table exactly match those in desired_columns
desired_columns_set = set(desired_columns)
table_columns_set = set(table.columns)

# Add missing columns to table
missing_columns = list(desired_columns_set - table_columns_set)

# Create a DataFrame with zeros for missing columns
missing_columns_df = pd.DataFrame(0, index=table.index, columns=missing_columns)

# Concatenate the missing columns DataFrame with the original table
table = pd.concat([table, missing_columns_df], axis=1)

# Remove extra columns from table
extra_columns = list(table_columns_set - desired_columns_set)
table = table.drop(columns=extra_columns)

# Reorder columns in table according to desired_columns
table = table[desired_columns]

#add empty phenotype column:
table["add_phenotype"] = np.nan

# Add genome_id column:
table.insert(0, 'genome_id', [genome_ids])

# Separate features (X) and target labels (y)
X_test = table.drop(['genome_id', 'add_phenotype'], axis=1)
y_test = LabelEncoder().fit_transform(table['add_phenotype'])

# Ensure that X_test is reshaped appropriately
X_test = X_test.values.reshape(X_test.shape[0], X_test.shape[1])

# Load the trained model
model = tf.keras.models.load_model(options.model)

# Make predictions on the test data
predictions = model.predict(X_test)
predicted_classes = (predictions > 0.5).astype(int)  # Adjust the threshold as needed

phenotype_mapping = {0: 'Resistant', 1: 'Susceptible'}
# Create a DataFrame with genome_id and predictions
result_df = pd.DataFrame({'genome_id': genome_ids, 'prediction': predicted_classes.flatten()})
result_df['predicted_phenotype'] = result_df['prediction'].map(phenotype_mapping)
# Write out the results file:
result_df.to_csv(options.output_file, index=False)
