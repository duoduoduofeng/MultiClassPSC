#!/bin/bash
import json
import sys
import os
import copy

from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score, f1_score
from sklearn.metrics import classification_report, confusion_matrix



def load_all_rs(rs_file):
	field_names = []
	all_sets = []

	with open(rs_file, 'r') as fin:
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
				cur_pro["real_distance"] = int(cur_pro["real_distance"])
				cur_pro["predict_class"] = int(cur_pro["predict_class"])

			all_sets.append(cur_pro)

	return all_sets


def calc_auc(all_sets, class_threshold, simi_threshold):
	true_labels = get_true_labels(all_sets, class_threshold)
	predicted_probs = get_probs(all_sets, class_threshold)

	rs = calc_metrics(true_labels, predicted_probs, simi_threshold)
	return rs


def get_true_labels(all_sets, class_threshold):
	true_labels = []

	for cur_pair in all_sets:
		distance = cur_pair["real_distance"]
		label = 1 if distance <= class_threshold else 0
		true_labels.append(label)

	return true_labels


def get_probs(all_sets, class_threshold):
	predicted_probs = []

	for cur_pair in all_sets:
		predict_distance = json.loads(cur_pair["predict_distance"])
		prob = 0
		for i in range(class_threshold + 1):
			prob += predict_distance[i]
		predicted_probs.append(prob)

	return predicted_probs


def calc_metrics(true_labels, predicted_probs, simi_threshold):
	predicted_labels = [1 if prob >= simi_threshold else 0 for prob in predicted_probs]

	accuracy = accuracy_score(true_labels, predicted_labels)
	precision = precision_score(true_labels, predicted_labels)
	recall = recall_score(true_labels, predicted_labels)
	auc = roc_auc_score(true_labels, predicted_probs)
	f1 = f1_score(true_labels, predicted_labels)

	rs = {}
	rs["accuracy"] = round(accuracy, 2)
	rs["precision"] = round(precision, 2)
	rs["recall"] = round(recall, 2)
	rs["auc"] = round(auc, 2)
	rs["f1"] = round(f1, 2)
	return rs


def calc_all_pr(all_sets):
	recall_dict = {}
	precision_dict = {}

	for cur_pair in all_sets:
		real_distance = cur_pair["real_distance"]
		discrete_predict_distance = cur_pair["predict_class"]

		isright = 0
		if real_distance == discrete_predict_distance:
			isright = 1

		# recall_dict statistics
		if real_distance not in recall_dict:
			recall_dict[real_distance] = {}

		if isright == 1:
			if "right" not in recall_dict[real_distance]:
				recall_dict[real_distance]["right"] = 0
			recall_dict[real_distance]["right"] += 1
		else:
			if "wrong" not in recall_dict[real_distance]:
				recall_dict[real_distance]["wrong"] = 0
			recall_dict[real_distance]["wrong"] += 1

		# precision_dict statistics
		if discrete_predict_distance not in precision_dict:
			precision_dict[discrete_predict_distance] = {}

		if isright == 1:
			if "right" not in precision_dict[discrete_predict_distance]:
				precision_dict[discrete_predict_distance]["right"] = 0
			precision_dict[discrete_predict_distance]["right"] += 1
		else:
			if "wrong" not in precision_dict[discrete_predict_distance]:
				precision_dict[discrete_predict_distance]["wrong"] = 0
			precision_dict[discrete_predict_distance]["wrong"] += 1

	whole_right = 0
	whole_wrong = 0
	for va in recall_dict:
		if "right" not in recall_dict[va]:
			recall_dict[va]["right"] = 0
		if "wrong" not in recall_dict[va]:
			recall_dict[va]["wrong"] = 0
		
		recall_dict[va]["recall"] = float(recall_dict[va]["right"]/(recall_dict[va]["right"] + recall_dict[va]["wrong"]))
		recall_dict[va]["recall"] = round(recall_dict[va]["recall"], 2)
		
		whole_right += recall_dict[va]["right"]
		whole_wrong += recall_dict[va]["wrong"]
	
	recall_dict["whole_recall"] = float(whole_right / (whole_right + whole_wrong))
	recall_dict["whole_recall"] = round(recall_dict["whole_recall"], 2)

	for va in precision_dict:
		if "right" not in precision_dict[va]:
			precision_dict[va]["right"] = 0
		if "wrong" not in precision_dict[va]:
			precision_dict[va]["wrong"] = 0
		
		precision_dict[va]["precision"] = float(precision_dict[va]["right"]/(precision_dict[va]["right"] + precision_dict[va]["wrong"]))
		precision_dict[va]["precision"] = round(precision_dict[va]["precision"], 2)

		whole_right += precision_dict[va]["right"]
		whole_wrong += precision_dict[va]["wrong"]
	
	precision_dict["whole_precision"] = float(whole_right / (whole_right + whole_wrong))
	precision_dict["whole_precision"] = round(precision_dict["whole_precision"], 2)
	
	return recall_dict, precision_dict


def calc_cr(all_sets):
	true_labels = []
	predicted_labels = []

	for cur_pair in all_sets:
		true_labels.append(cur_pair["real_distance"])
		predicted_labels.append(cur_pair["predict_class"])

	# print("\nClassification Report:")
	# print(classification_report(true_labels, predicted_labels))

	cr = classification_report(true_labels, predicted_labels, output_dict=True)
	return cr



def display_avg(rs_dir, algo = "PDM", scenario = "bothin", class_threshold = 1):
	metrics = {}
	
	for root, dirs, files in os.walk(rs_dir):
		for thefile in files:
			all_sets = load_all_rs(os.path.join(root, thefile))

			cur_metrics = {}
			max_f1 = -1.0
			for simi_threshold in [i / 10 for i in range(1, 10)]:
				rs = calc_auc(all_sets, class_threshold, simi_threshold)

				f1 = rs["f1"]
				if f1 <= max_f1:
					continue

				# Choose the one with maximum f1 score.
				max_f1 = f1
				cur_metrics["f1"] = f1
				for met in ["auc", "accuracy", "precision", "recall"]:
					cur_metrics[met] = rs[met]
				cur_metrics["simi_threshold"] = simi_threshold

			for met in ["auc", "accuracy", "precision", "recall", "f1", "simi_threshold"]:
				if met not in metrics:
					metrics[met] = []
				metrics[met].append(cur_metrics[met])


	avg_metrics = {}
	for ele in metrics:
		if len(metrics[ele]) == 0:
			avg_metrics[ele] = 0
		avg_metrics[ele] = sum(metrics[ele]) / len(metrics[ele])
		avg_metrics[ele] = round(avg_metrics[ele], 2)

	auc = avg_metrics["auc"]
	# simi_threshold = avg_metrics["simi_threshold"]
	accuracy = avg_metrics["accuracy"]
	precision = avg_metrics["precision"]
	recall = avg_metrics["recall"]
	f1 = avg_metrics["f1"]
	print(f" & {algo} & {auc} & {accuracy} & {precision} & {recall} & {f1} \\\\")

	return avg_metrics



def display_single(rs_file, algo = "PDM", class_threshold = 1):
	all_sets = load_all_rs(rs_file)

	if algo == "PDM":
		recall_dict, precision_dict = calc_all_pr(all_sets)
		print("Precision and Recall for all categories:")
		print(json.dumps(recall_dict))
		print(json.dumps(precision_dict))
	
	print("\n\nMetrics of binary module:")
	for simi_threshold in [i / 10 for i in range(1, 10)]:
		rs = calc_auc(all_sets, class_threshold, simi_threshold)
		print(json.dumps(rs))

		auc = rs["auc"]
		accuracy = rs["accuracy"]
		precision = rs["precision"]
		recall = rs["recall"]
		f1 = rs["f1"]
		print(f" & {algo} & {auc} & {simi_threshold} & {accuracy} & {precision} & {recall} & {f1} \\\\")


def display_single_cr(rs_file):
	all_sets = load_all_rs(rs_file)
	cr = calc_cr(all_sets)

	# print(json.dumps(cr))
	return cr



def display_avg_cr(rs_dir):
	metric_dicts = []
	for root, dirs, files in os.walk(rs_dir):
		for thefile in files:
			rs_file = os.path.join(root, thefile)
			print(f"Dealing {rs_file}.")
			thecr = display_single_cr(rs_file)
			metric_dicts.append(copy.deepcopy(thecr))


	avg_metrics = {}
	for cl in ["1", "2", "3", "4", "weighted avg"]:
		for met in ["precision", "recall", "f1-score"]:
			vs = []
			for thedict in metric_dicts:
				vs.append(thedict[cl][met])
			if cl not in avg_metrics:
				avg_metrics[cl] = {}
			avg_metrics[cl][met] = round(sum(vs) / len(vs), 2)

	# accuracy
	vs = []
	for thedict in metric_dicts:
		vs.append(thedict["accuracy"])
	avg_metrics["accuracy"] = round(sum(vs) / len(vs), 2)

	# print(json.dumps(avg_metrics))
	return avg_metrics


if __name__ == "__main__":
	# rs_file = sys.argv[1]
	
	subset_name = "cl_1000001"
	base_dir = "/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/evaluation_result"
	
	scenario = "onlyone"
	validation_size = 1000
	# rs_file = f"{base_dir}/{subset_name}/result/{validation_size}/pdms/{scenario}/8.txt"
	
	###### Binary Classification
	# for validation_size in [200, 400, 600, 800, 1000, 2000, 3000, 10000]:
	# print("\n\n***************************************************")
	# for scenario in ["bothin", "onlyone", "neither"]:
	# 	rs_dir = f"{base_dir}/{subset_name}/result/{validation_size}/pdms/{scenario}"
	# 	print(f"\n==== {validation_size} ==== {scenario} ====")
	# 	print(f"{rs_dir}")
	# 	display_avg(rs_dir, algo = "PDM", scenario = scenario, class_threshold=2)
	# print("\n")


	###### Multi-Classification
	# single
	# display_single_cr(rs_file)

	# batch
	rs_dir = f"{base_dir}/{subset_name}/result/{validation_size}/pdms/{scenario}"
	display_avg_cr(rs_dir)



