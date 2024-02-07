# -Align-COVID-19-Genomes-and-Analyze-Mutations
#Analyze COVID-19 genomes for mutations, aligning with NC_045512. Perform statistical analysis, plot mutation frequency, SNP types, and Indels. Use negative binomial distribution for frequent mutations. Calculate mutation matrix and MST with GATEkeeper. Validate MST order against sequencing time
def mainGenomeToFile(genomes, x):
    file1 = open("f3.fasta", "w")
    file1.write(">" + genomes[x].id + "\n")
    file1.write(str(genomes[x].seq))
    file1.close()
 
def recToFile(genomes, i1):
    file1 = open("f4.fasta", "w")
    file1.write(">" + genomes[i1].id + "\n")
    file1.write(str(genomes[i1].seq))
    file1.close()
 
import os
 
def GATEkeeperAlign():
    command = "bin/GATEkeeper -r f3.fasta -q f4.fasta -o output2"
    os.system(command)
 
def countMutations():
    fileB = open("output2.vcf", "r")
    mutationsCount = 0
    for line in fileB:
        if line.startswith("#"):
            continue
        else:
            mutationsCount += 1
    fileB.close()
    return mutationsCount
 
def countIndel(info):
    deletions = info[info['INFO'].str.endswith('DELETE')]
    insertions = info[info['INFO'].str.endswith('INSERT')]
    delCount = len(deletions)
    insCount = len(insertions)
    IndelCount = delCount + insCount
    return IndelCount, insCount, delCount
 
def countSNP(info):
    substitutions = info[info['INFO'].str.endswith('SUBSTITUTE')]
    SNPcount = len(substitutions)
    return SNPcount
 
def SNPtype(info):
    substitutions = info[info['INFO'].str.endswith('SUBSTITUTE')]
    SNPtypes = {}
    for nt in "AGCT":
        subset = substitutions[substitutions['REF'] == nt]
        for (v1, v2) in subset['ALT'].value_counts().items():
            SNPtypes[f"{nt}-{v1}"] = v2
    return SNPtypes
 
import pandas as pd
 
def getRows():
    file1 = open("output2.vcf", "r")
    vcfRows = []
    for line in file1:
        if line.startswith("#"):
            continue
        else:
            vcfRows.append(line.strip().split('\t'))
    file1.close()
    return vcfRows
 
import Bio
from Bio import SeqIO
from Bio import SeqRecord
import numpy as np
 
genomes1000 = []
count = 0
adjMatrix = np.zeros((1, 999))
 
for record in SeqIO.parse("Selected_Unique_COVID19_Genomes_Asia.fasta", "fasta"):
    if record.id == "NC_045512":
        genomes1000.append(record)
        count += 1
 
for record in SeqIO.parse("Selected_Unique_COVID19_Genomes_Asia.fasta", "fasta"):
    genomes1000.append(record)
    count += 1
    if count == 1000:
        break
 
genomeNames1000 = []
for item in genomes1000:
    genomeNames1000.append(item.id)
 
i = genomeNames1000.index("NC_045512.2")
 
mainGenomeToFile(genomes1000, i)
 
rows = []
totalRows = []
 
for val in range(1, 999):
    recToFile(genomes1000, val)
    GATEkeeperAlign()
    rows.extend(getRows())
    adjMatrix[0, val] = countMutations()
 
vcfInfo = pd.DataFrame(totalRows, columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
 
indelNumber, insNumber, delNumber = countIndel(vcfInfo)
countSNP = countSNP(vcfInfo)
types = SNPtype(vcfInfo)
 
vcfInfo = pd.read_csv("AllVCF999.csv")
 
indelNumber, insNumber, delNumber = countIndel(vcfInfo)
countSNP = countSNP(vcfInfo)
types = SNPtype(vcfInfo)
 
import seaborn as sns
import matplotlib.pyplot as plt
 
counts = [countSNP, delNumber, insNumber]
sns.barplot(x=['Substitutions', 'Deletions', 'Insertions'], y=counts)
for i in range(3):
    plt.text(i, counts[i] + 10, str(counts[i]), ha='center')
 
plt.title("Number of each mutation type")
 
import seaborn as sns
 
sns.barplot(x=list(types.keys()), y=list(types.values()))
for i, (key, val) in enumerate(types.items()):
    plt.text(i, val + 30, val, ha='center')
 
plt.title("SNP Type Frequency")
 
vcfInfo.to_csv("AllVCF999.csv")
 
import pandas as pd
import matplotlib.pyplot as plt
 
info = vcfInfo
 
deletions = info[info['INFO'].str.endswith('DELETE')]
insertions = info[info['INFO'].str.endswith('INSERT')]
 
indels = pd.concat([insertions, deletions])
indels['REF'].str.len().value_counts().plot(kind='bar', color="green")
plt.title("Analysis of Indel Frequency")
plt.xlabel("Number of Nucleotides Inserted/Deleted")
plt.ylabel("Frequency")
 
import matplotlib.pyplot as plt
 
posCount = vcfInfo.groupby(['POS'], as_index=False).count()
 
positions = posCount['POS'].values.tolist()
frequencies = posCount["ID"].values.tolist()
 
plt.scatter(positions, frequencies, s=10, color='blue')
plt.xlabel("Position on Chromosome")
plt.ylabel("Frequency of Mutations")
plt.title("Analysis of Mutation Frequency at Each Position")
plt.show()
 
import pandas as pd
 
df = posCount
 
value_to_find = 995
 
result = df[df['ID'] == value_to_find]
 
if not result.empty:
    print(f"The value '{value_to_find}' is present in the DataFrame at row:")
    print(result)
else:
    print(f"The value '{value_to_find}' is not found in the DataFrame.")
 
posCount
n = 0
 
posCount
 
df = posCount
 
for row in posCount:
    if n == posCount['POS'].iloc[n]:
        continue
    else:
        new_row = pd.DataFrame({'POS': [n], 'Unnamed: 0': [0], '#CHROM': [0],
                                'ID': [0], 'REF': [0], 'ALT': [0], 'QUAL': [0],
                                'FILTER': [0], 'INFO': [0]})
    insert_index = n
    df = pd.concat([df.iloc[:insert_index], new_row, df.iloc[insert_index:]]).reset_index(drop=True)
 
    n += 1
    if n == 29841:
        break
 
df
frequencies = df["ID"].values.tolist()
 
from scipy.optimize import minimize
from scipy.stats import nbinom
import numpy as np
 
mutation_counts = frequencies
 
def neg_log_likelihood(params, data):
    r, p = params
    return -np.sum(nbinom.logpmf(data, r, p))
 
mean_count = np.mean(mutation_counts)
var_count = np.var(mutation_counts)
 
initial_guess = [mean_count, 0.5]
result = minimize(neg_log_likelihood, initial_guess, args=(mutation_counts,), method='Nelder-Mead')
 
estimated_r, estimated_p = result.x
print("Estimated parameters (r, p):", estimated_r, estimated_p)
 
from scipy.stats import nbinom
import matplotlib.pyplot as plt
import numpy as np
 
mutation_frequencies = frequencies
 
p_values = [nbinom.sf(count, estimated_r, estimated_p) for count in mutation_frequencies]
 
significance_threshold = 0.0001
 
significant_positions = [i for i, p_value in enumerate(p_values) if p_value < significance_threshold]
 
print("Significant positions:", significant_positions)
 
import matplotlib.pyplot as plt
from scipy.stats import nbinom
import numpy as np
 
mutation_counts = frequencies
estimated_r = estimated_r
estimated_p = estimated_p
 
x_values = np.arange(0, np.max(mutation_counts) + 1)
 
pmf = nbinom.pmf(x_values, estimated_r, estimated_p)
 
plt.bar(x_values, pmf, color='skyblue', label='Fitted Distribution')
plt.xlabel('Number of Mutations')
plt.ylabel('Probability')
plt.title('Fitted Negative Binomial Distribution')
plt.ylim(0, 0.1)
plt.legend()
plt.show()
 
import scipy as sp
from scipy import sparse
 
def generateMST(data):
    adjMatrixB = sp.sparse.csr_matrix(data)
    mst = sp.sparse.csgraph.minimum_spanning_tree(adjMatrixB, overwrite=False)
    return mst
 
import csv
import numpy as np
import pandas as pd
 
adjMatrix = []
adj1000 = pd.read_csv("AdjacencyMatrix1000.csv", index_col=0)
adjMatrix = adj1000.to_numpy()
 
mst = generateMST(adjMatrix)
 
mst1 = mst.toarray().astype(int)
 
np.savetxt("MST1000.csv", mst1, delimiter=",", fmt='%a')
 
import numpy as np
from scipy.sparse import csgraph
import networkx as nx
import matplotlib.pyplot as plt
 
graph_matrix = mst1
 
mst2 = csgraph.minimum_spanning_tree(graph_matrix)
 
G = nx.from_scipy_sparse_matrix(mst2)
 
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, node_size=100, node_color='skyblue', font_weight='light', font_size=5)
plt.title('Minimum Spanning Tree For 1000 Adjacency Matrix')
plt.savefig('MST1000Visual.png', format='png', dpi=600)
plt.show()

