import numpy as np
import pandas as pd
from collections import defaultdict
import scipy.stats as stat
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression

def reactome_genes():
	output = defaultdict(list)
	output_list = []
	f = open('./data/Sour_Data/c2.all.v7.2.symbols.gmt','r')
	lines = f.readlines()
	for line in lines:
		line = line.strip().split('\t')
		if 'REACTOME' in line[0]:
			reactome = line[0]
			output_list.append(reactome)
			for i in range(2, len(line)):
				gene = line[i]
				output[reactome].append(gene)
	f.close()
	return output


if __name__ == "__main__":
	all_genes = pd.read_csv("./data/Sour_Data/geneid_info.txt", sep="\t")
	reactome = reactome_genes()
	
	propagated_network_df = pd.read_csv("./data/Result/propagation_score.txt", sep="\t")
	propagated_network_df = propagated_network_df.sort_values("Propagate_Score", ascending= False)
    # 选择前150个基因来剪裁蛋白质相互作用网络
	selected_genes = propagated_network_df.iloc[1:100]["Gene_Symbol"].to_list()

	q_cutoff = 0.005
	i = 0
	temp_hypergeom = defaultdict(list)
	M = len(all_genes) # M代表所有基因的个数
	N = len(selected_genes) # N代表筛选的100个基因
	pvals, qvals, overlap_counts, pw_counts = [], [], [], [] # 保存pval，qval，筛选基因与通路基因重叠个数，通路基因个数
	for pw_name in list(reactome):
		pw_genes = list(set(reactome[pw_name]) & set(all_genes['Gene_Symbol'])) # 通路基因名称
		n = len(pw_genes) # n代表通路基因个数
		k = len(set(pw_genes) & set(selected_genes)) # k代表通路基因与筛选的100个基因的重叠基因数

		p = stat.hypergeom.sf(k, M, n, N)  # 超几何分布
		pvals.append(p)
		overlap_counts.append(k)
		pw_counts.append(n)
		temp_hypergeom['pw_name'].append(pw_name)
		temp_hypergeom['p'].append(p)
		if (i % 100 == 1):
			# print(p)
			# print(k)
			# print(pw_genes)
			pass
		i = i + 1
	reject_, qvals, alphacSidak, alphacBonf = multipletests(pvals)
	temp_hypergeom['qval'] = qvals
	temp_hypergeom = pd.DataFrame(temp_hypergeom)
	pw_features = temp_hypergeom.loc[temp_hypergeom['qval'] < q_cutoff, :]['pw_name'].tolist()
	pd.DataFrame(pw_features).to_csv("./data/Result/pw_feature.txt", index=False, header=False)
	print(pw_features)
