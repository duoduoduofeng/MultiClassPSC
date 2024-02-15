from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.Align import substitution_matrices
from sklearn.metrics import precision_score, recall_score, roc_auc_score

import os
import json


def calc_whole_blast_sim(validate_file, blast_sim_file):
	field_names = []

	# 保存序列到文件
	query_seq_file = "query.fasta"
	subject_seq_file = "subject.fasta"
	output_file = "blast_output.xml"

	with open(validate_file, 'r') as fin, open(blast_sim_file, 'w') as fout:
		row_num = 0
		for line in fin:
			row_num += 1
			cur_pro = {}

			parts = line.strip().split("\t")

			if row_num == 1:
				field_names = parts
				fout.write(f"{line.strip()}\tblast_simi\tnw_simi\tsw_simi\n")
				continue
			else:
				for i in range(0, len(parts)):
					cur_pro[field_names[i]] = parts[i]
				cur_pro["distance"] = int(cur_pro["distance"])

			query_seq = cur_pro["protein1_seq"]
			subject_seq = cur_pro["protein2_seq"]

			with open(query_seq_file, "w") as f:
				f.write(">query\n" + query_seq)
			with open(subject_seq_file, "w") as f:
				f.write(">subject\n" + subject_seq)
			
			calc_blast_sim(query_seq_file, subject_seq_file, output_file)
			blast_sim = parse_record(output_file, row_num)
			global_similarity, local_similarity = \
				calc_nw_sw(query_seq, subject_seq)

			cur_pro["blast_simi"] = blast_sim
			cur_pro["nw_simi"] = global_similarity
			cur_pro["sw_simi"] = local_similarity

			fout.write(f"{line.strip()}\t{blast_sim}\t{global_similarity}\t{local_similarity}\n")

	# clear
	os.remove(output_file)
	os.remove(query_seq_file)
	os.remove(subject_seq_file)



def load_simis(blast_sim_file):
	field_names = []
	all_sets = []

	with open(blast_sim_file, 'r') as fin:
		row_num = 0
		for line in fin:
			row_num += 1
			cur_pro = {}

			parts = line.strip().split("\t")

			if row_num == 1:
				field_names = parts
				continue
			else:
				for i in range(0, len(parts)):
					cur_pro[field_names[i]] = parts[i]
				cur_pro["distance"] = int(cur_pro["distance"])
				cur_pro["blast_simi"] = float(cur_pro["blast_simi"])
				cur_pro["nw_simi"] = float(cur_pro["nw_simi"])
				cur_pro["sw_simi"] = float(cur_pro["sw_simi"])

			all_sets.append(cur_pro)

	# Calculate similarity.
	rs = {}
	# class_threshold = 1
	# simi_threshold = 0.4
	true_labels = get_true_labels(all_sets, class_threshold=1)
	rs["blast"] = blast_metrics(all_sets, true_labels, simi_threshold=0.5)
	rs["nw"] = nw_metrics(all_sets, true_labels, simi_threshold=0.4)
	rs["sw"] = sw_metrics(all_sets, true_labels, simi_threshold=0.1)

	print("Metrics of different benchmark.")
	print(json.dumps(rs))


def calc_blast_sim(query_seq_file, subject_seq_file, output_file):
	# 设置BLAST命令行
	blastp_cline = NcbiblastpCommandline(query=query_seq_file, 
		subject=subject_seq_file, outfmt=5, out=output_file)

	# 执行BLAST
	stdout, stderr = blastp_cline()


def parse_record(output_file, flag):
	with open(output_file, 'r') as fin:
		blast_record = NCBIXML.read(fin)

		# 存储所有比对的身份值和比对长度
		identities = []
		align_lengths = []

		# 获取所有比对的身份值和比对长度
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				alignment_length = hsp.align_length
				identity = (hsp.identities / alignment_length)
				identities.append(identity)
				align_lengths.append(alignment_length)

		# 计算加权平均的 per identity
		if identities:
			total_identity = sum(identity * align_length for identity, align_length in zip(identities, align_lengths))
			total_align_length = sum(align_lengths)
			average_per_identity = total_identity / total_align_length
			average_per_identity = round(average_per_identity, 2)
			# print("Average Per Identity:", average_per_identity)
			return average_per_identity
		else:
			# print(f"Line {flag}, No alignment results found")
			return 0


def calc_nw_sw(seq1, seq2):
	# 获取BLOSUM62打分矩阵
	matrix = substitution_matrices.load('BLOSUM62')

	# 全局对齐 - Needleman-Wunsch算法
	global_alignment = pairwise2.align.globalds(seq1, seq2, matrix, -10, -0.5)[0]  # gap开启和扩展惩罚

	# 局部对齐 - Smith-Waterman算法
	local_alignment = pairwise2.align.localds(seq1, seq2, matrix, -10, -0.5)[0]

	# 计算相似度值
	def calculate_similarity(alignment, length):
	    score, max_score = alignment.score, length * matrix[('A', 'A')]  # 假设最大得分为序列长度乘以最高单个匹配得分
	    similarity = score / max_score
	    return similarity

	# 计算全局和局部对齐的相似度
	global_similarity = calculate_similarity(global_alignment, min(len(seq1), len(seq2)))
	local_similarity = calculate_similarity(local_alignment, min(len(seq1), len(seq2)))
	global_similarity = round(global_similarity, 2)
	local_similarity = round(local_similarity, 2)

	return global_similarity, local_similarity


def blast_metrics(all_sets, true_labels, simi_threshold):
	predicted_probs = get_probs(all_sets, "blast_simi")
	return cal_metrics(true_labels, predicted_probs, simi_threshold)


def nw_metrics(all_sets, true_labels, simi_threshold):
	predicted_probs = get_probs(all_sets, "nw_simi")
	normalized_predicted_probs = max_norm(predicted_probs)
	return cal_metrics(true_labels, normalized_predicted_probs, simi_threshold)


def sw_metrics(all_sets, true_labels, simi_threshold):
	predicted_probs = get_probs(all_sets, "sw_simi")
	normalized_predicted_probs = max_norm(predicted_probs)
	return cal_metrics(true_labels, normalized_predicted_probs, simi_threshold)


def max_norm(data):
	# 找到列表中的最大值和最小值
	min_val = min(data)
	max_val = max(data)

	# 对列表进行最大-最小归一化处理
	normalized_data = [(x - min_val) / (max_val - min_val) for x in data]

	return normalized_data


def get_true_labels(all_sets, class_threshold):
	true_labels = []

	for cur_pair in all_sets:
		distance = cur_pair["distance"]
		label = 1 if distance <= class_threshold else 0
		true_labels.append(label)

	return true_labels


def get_probs(all_sets, sim_filed):
	predicted_probs = []

	for cur_pair in all_sets:
		predicted_probs.append(cur_pair[sim_filed])

	return predicted_probs



def cal_metrics(true_labels, predicted_probs, simi_threshold):
	predicted_labels = [1 if prob >= simi_threshold else 0 for prob in predicted_probs]

	precision = precision_score(true_labels, predicted_labels)
	recall = recall_score(true_labels, predicted_labels)
	auc = roc_auc_score(true_labels, predicted_probs)

	rs = {}
	rs["precision"] = round(precision, 2)
	rs["recall"] = round(recall, 2)
	rs["auc"] = round(auc, 2)
	return rs


if __name__ == "__main__":
	# validate_file = "/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/datasets/trial/data/sample_proteins_dataset.validate.txt"
	blast_sim_file = f"/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/src/evaluates/test.blast_sim"
	
	validate_file = "/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/model_bak/MultiClassPSC/generated_data/datasets/try_cl_1000001.5class/data/20240130_205512.best_data/sample_proteins_dataset.validate.txt.best"

	# calc_whole_blast_sim(validate_file, blast_sim_file)
	load_simis(blast_sim_file)


