import pandas as pd

R_regulate_network = pd.read_csv("./data/Hand_Data/R_regulate_network.txt", sep=",")
NR_regulate_network = pd.read_csv("./data/Hand_Data/NR_regulate_network.txt", sep=",")

gene_list = list(set(R_regulate_network['Gene1'].tolist()) | set(R_regulate_network["Gene2"].tolist()))
expr_genes1 = pd.read_csv("./data/Sour_Data/Kim/Kim_log2tpm.txt", sep="\t")
expr_genes2 = pd.read_csv("./data/Sour_Data/GSE91061/GSE91061_log2tpm.txt", sep="\t")
expr_genes3 = pd.read_csv("./data/Sour_Data/GSE115821/GSE115821_log2tpm.txt", sep="\t")

selected_genes = list(set(gene_list) & set(expr_genes1.loc[:,'gene_id'].tolist()) & set(expr_genes2.loc[:,'gene_id'].tolist())
                      & set(expr_genes3.loc[:,'gene_id'].tolist()))
# 给基因编码
gene_dict = {}
for i in range(len(selected_genes)):
    gene_dict[selected_genes[i]] = i

# 生成edge_index数组 [2, n]，n代表边的个数
edge_index = [[],[]]
for i in range(len(R_regulate_network)):
    gene1 = R_regulate_network.iloc[i, 0]
    gene2 = R_regulate_network.iloc[i, 1]
    corVal = R_regulate_network.iloc[i,2]
    if gene_dict.__contains__(gene1) and gene_dict.__contains__(gene2) and abs(corVal) > 0.145:
        edge_index[0].append(gene_dict[gene1])
        edge_index[1].append(gene_dict[gene2])
print(len(edge_index[0]))