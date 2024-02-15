import sys
import json
import random
import copy
import math
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
from collections import Counter


"""
Count the frequency of each {pdb_id, chain_id} pair.
"""
def uniq_proteins_count(proteins_file):
	# Dict runs much faster than list
	# {pdb_id: {chain_id, count}}
	uniq_pros = {}

	with open(proteins_file["filename"], 'r') as f2:
		row_num = 0
		for line in f2:
			row_num += 1
			if row_num < proteins_file["start_row_num"]:
				continue

			parts = line.strip().split(' ')
			if len(parts) == proteins_file["fields_count"]:
				pdb_id = parts[1]
				chain_id = parts[2]
				if pdb_id not in uniq_pros:
					new_chains = {}
					new_chains[chain_id] = 1
					uniq_pros[pdb_id] = new_chains
				else:
					cur_chain = uniq_pros[pdb_id]
					if chain_id not in cur_chain:
						cur_chain[chain_id] = 1
					else:
						cur_chain[chain_id] = cur_chain[chain_id] + 1

	print(f"There are {len(uniq_pros)} uniq pdb_id in total.")

	return uniq_pros


def to_count(uniq_pros):
	total_uniq_pdb_chain = 0
	for pdb_id in uniq_pros:
		total_uniq_pdb_chain += len(uniq_pros[pdb_id])

	print(f"Uniq pdb_id, chain_id pairs count: {total_uniq_pdb_chain}")



"""
Reserve the {pdb_id, chain_id} pairs which occurs exactly twice.
PS, only family representative domains are keeped. 
DOMID REPRESENTED-PDBID REPRESENTED-PDBCHAIN FA-DOMID FA-PDBID FA-PDBREG FA-UNIID FA-UNIREG SF-DOMID SF-PDBID SF-PDBREG SF-UNIID SF-UNIREG SCOPCLA
8000061 2DT5 B 8000061 2DT5 A:2-77 Q5SHS3 2-77 8001519 2DT5 A:4-63 Q5SHS3 4-63 TP=1,CL=1000003,CF=2000145,SF=3000034,FA=4000057
"""
def filter_multi_domains_protein(uniq_pros, uniq_pros_previous_file, merged_info_file):
	# generate an assist file, to help locate the sequence filename.
	row_flags = {}
	with open(uniq_pros_previous_file, 'r') as fin:
		row_num = 0
		for line in fin:
			pdb_id = line.split("\t")[0]
			row_num += 1
			row_flags[pdb_id] = row_num

	reserved_proteins = {}
	with open(merged_info_file, 'r') as fin:
		row_num = 0
		sf_row_count = 0
		for line in fin:
			row_num += 1
			if row_num < 2:
				continue

			parts = line.strip().split(" ")
			if len(parts) == 14:
				rep_id = parts[0]
				pdb_id = parts[1]
				chain_id = parts[2]
				lineage_info = parts[-1]
				# if uniq_pros[pdb_id][chain_id] > -1:
				# Here filters the proteins which contains not 2 domains.
				if uniq_pros[pdb_id][chain_id] == 2:
					thekey = f"{pdb_id}_{chain_id}"
					info = {}
					info["rep_id"] = rep_id
					info["sequence_file_row_num"] = row_flags[pdb_id]
					# TP=1,CL=1000002,CF=2000148,SF=3000038,FA=4000119
					for thestr in lineage_info.split(","):
						nodename, nodevalue = thestr.split("=")
						info[nodename] = nodevalue
					reserved_proteins[thekey] = info

			elif len(parts) == 3:
				# print(f"**** Line {row_num} is not for family domain but super family's.")
				sf_row_count += 1
				continue
			else:
				print(f"Invalid input, line {row_num} in {merged_info_file}.")
				continue

	print(f"There are {len(reserved_proteins)} proteins reserved. \nAnd there are {sf_row_count} super family lines.\nPrint some samples to check the reserved proteins.")
	# an_count = 0
	# for thekey in reserved_proteins:
	# 	an_count += 1
	# 	if an_count%10000 == 1:
	# 		print(f"{thekey}\t{json.dumps(reserved_proteins[thekey])}")

	return reserved_proteins


"""
Stat the lineage distribution of the reserved proteins.
"""
def stat_lineage(reserved_proteins):
	stats = {
		"TP": {}, 
		"CL": {},
		"CF": {},
		"SF": {}, 
		"FA": {}, 
		"rep_id": {}
	}
	for thekey in reserved_proteins:
		info = reserved_proteins[thekey]
		if info["CL"] != "1000001":
			continue
		for node in ["rep_id", "TP", "CL", "CF", "SF", "FA"]:
			if info[node] not in stats[node]:
				stats[node][info[node]] = 1
			else:
				stats[node][info[node]] = stats[node][info[node]] + 1
	
	# print(json.dumps(stats))

	for nodename in ["CF", "SF", "FA", "rep_id"]:
		# plot_dis(stats, nodename)
		# plot_dis(stats, nodename, 100)
		# plot_dis(stats, nodename, 1000)
		plot_scatter(stats, nodename)

	# return stats


def plot_dis(stats, nodename, th = 0):
	# Data from the list	
	data = list(stats[nodename].values())
	# if filter
	if th > 0:
		data = [x for x in data if x <= th]
	# print(len(data))
	# print(data)
	# print()

	# thebins = [0, 20, 50, 100, 200, 500, 1000, 2000, 5000, 7000]

	# Creating bar chart
	plt.clf()
	plt.hist(data, bins=50, color='skyblue', edgecolor='black')
	# plt.xticks(thebins)
	plt.xlabel('Value')
	plt.ylabel('Frequency')
	plt.title(f"Distribution of protein count of different {nodename} in CL=1000001.")

	# Saving the bar chart as an image
	bar_chart_name = f"../../generated_data/statistics/{nodename}.cl_1000001.bar_chart.png"
	if th > 0:
		bar_chart_name = f"../../generated_data/statistics/{nodename}.cl_1000001.lt_{th}.bar_chart.png"
	plt.savefig(bar_chart_name)

	# Calculating statistical measures
	mean_value = np.mean(data)
	median_value = np.median(data)
	mode_value = Counter(data).most_common(1)[0][0]
	std_dev = np.std(data)
	variance = np.var(data)
	data_range = np.max(data) - np.min(data)
	minv, q1, q3, maxv = np.percentile(data, [0, 25, 75, 100])

	# Displaying statistical measures
	stat_file = f"../../generated_data/statistics/{nodename}.cl_1000001.stats"
	if th > 0:
		stat_file = f"../../generated_data/statistics/{nodename}.cl_1000001.lt_{th}.stats"
	with open(stat_file, 'w') as fout:
		fout.write(f"Mean: {mean_value}\n")
		fout.write(f"Median: {median_value}\n")
		fout.write(f"Mode: {mode_value}\n")
		fout.write(f"Standard Deviation: {std_dev}\n")
		fout.write(f"Variance: {variance}\n")
		fout.write(f"Range: {data_range}\n")
		fout.write(f"Minimum: {minv}\n")
		fout.write(f"1st Quartile (Q1): {q1}\n")
		fout.write(f"3rd Quartile (Q3): {q3}\n")
		fout.write(f"Maximum: {maxv}\n")


def plot_scatter(stats, nodename, th = 0):
	# Data from the list	
	data = list(stats[nodename].values())
	# if filter
	if th > 0:
		data = [x for x in data if x <= th]

	# 使用Counter类计数
	counter = Counter(data)

	# 将计数结果分解成x和y坐标
	x_values = sorted(counter.keys())
	y_values = [counter[x] for x in x_values]

	# 绘制折线图
	plt.plot(x_values, y_values, linestyle='-')

	# 添加标签和标题
	plt.xlabel('Value')
	plt.ylabel('Count')
	plt.title('Line Plot of Data Counts')

	# 显示图形
	# Saving the bar chart as an image
	sacatter_chart_name = f"../../generated_data/statistics/{nodename}.cl_1000001.scatter.png"
	if th > 0:
		sacatter_chart_name = f"../../generated_data/statistics/{nodename}.cl_1000001.lt_{th}.scatter.png"
	plt.savefig(sacatter_chart_name)



if __name__ == "__main__":
	meta_info = [
	    {
	        "filename": "../../data/scop-represented-structures-latest.txt",
	        "start_row_num": 7,
	        "fields_count": 3,
	        "head": "DOMID REPRESENTED-PDBID REPRESENTED-PDBCHAIN",
	        "example": "8000061 2DT5 B"
	    }
	]

	### Step 1, stat all the proteins
	proteins_file = meta_info[0]
	uniq_pros = uniq_proteins_count(proteins_file)

	# to_count(uniq_pros)

	### Step 2, 
	uniq_pros_previous_file = "../../generated_data/whole_pdbs/uniq_pdb_list.txt"
	merged_info_file = "../../generated_data/scop-protein-whole-info.txt"
	reserved_proteins = filter_multi_domains_protein(uniq_pros, uniq_pros_previous_file, merged_info_file)

	### Step 3, stat the reserved proteins lineage.
	stats = stat_lineage(reserved_proteins)
