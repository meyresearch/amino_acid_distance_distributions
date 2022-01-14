import pandas as pd
import numpy as np

DF = pd.read_csv('../data/duplicates.csv')

#DF_T = DF.transpose()

string_list = DF['DOWNLOADED'].astype(str).to_numpy()

down_list = [x.replace('.pdb', '').replace('.gz', '') for x in string_list]
all_list = DF['ALL'].astype(str).tolist()

duplicates = []

for i in range(len(all_list)):
    print(f'We are at entry: {i}')
    duplicates.append(down_list[i])
    duplicates.append(all_list[i])


duplicate_dict = {'duplicates':duplicates}
print(f'Length of duplicates list is {len(duplicates)}')

df = pd.DataFrame()
df = df.from_dict(duplicate_dict)

no_duplicates = df.drop_duplicates(inplace=False)

no_dupes_list = no_duplicates['duplicates'].to_numpy()
print(f'No duplicates list length is {len(no_dupes_list)}')

np.savetxt('no_duplicates.csv',no_dupes_list,delimiter=',')
#print(no_duplicates)