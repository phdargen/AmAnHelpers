import uproot
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

#file_path = "files/psi_data_bkgPdf.root"
file_path = "files/psi_mc_pidCorr_rw.root"
tree_name = "DecayTree"

with uproot.open(file_path) as file:
    tree = file[tree_name]
    array = tree.arrays(library="np")

    # Then create a DataFrame from the NumPy array
    df = pd.DataFrame(array)

print(df.head())
branch_names = df.columns.tolist()
print("Branch Names:")
print(branch_names)

# ----

# Cell 4: Visualize the data using Seaborn

# Replace with the column name from your DataFrame that you'd like to visualize
#column_to_visualize = "column_name"

# Histogram
#sns.histplot(df[column_to_visualize])
#plt.title(f"Histogram of {column_to_visualize}")
#plt.show()

variables = ["B_DTF_MM", "psi_PV_MM", "DTF_mKpipi", "DTF_mKpi", "DTF_mpsipipi", "DTF_mpsipi"]
sns.set(rc={'figure.figsize':(10, 10)})  # You can adjust the size as needed
g = sns.pairplot(df, vars=variables, height=2)  # The height parameter controls the size of the subplots
plt.savefig("data.png", dpi=300, bbox_inches='tight')
