import json
import matplotlib.pyplot as plt
import numpy as np

import benchmarks
import evaluate


def obtain_metric_values(base_dir, subset_name):
	all_mets = {}
	
	for validation_size in [200, 400, 600, 800, 1000, 2000, 3000]:
		for scenario in ["bothin", "onlyone", "neither"]:
			print(f"\n========== {validation_size} ========== {scenario} ===========")
			blast_sim_dir = f"{base_dir}/{subset_name}/result/{validation_size}/benchmarks/{scenario}"
			
			cur_metrics = {}

			cur_metrics["nw"], cur_metrics["sw"], cur_metrics["blast"] = \
				benchmarks.display_avg(blast_sim_dir, class_threshold=2)

			rs_dir = f"{base_dir}/{subset_name}/result/{validation_size}/pdms/{scenario}"
			cur_metrics["pdm"] = evaluate.display_avg(rs_dir, algo = "PDM", scenario = scenario, class_threshold=1)

			for model in ["pdm", "nw", "sw", "blast"]:
				if model not in all_mets:
					all_mets[model] = {}
				if validation_size not in all_mets[model]:
					all_mets[model][validation_size] = {}
				if scenario not in all_mets[model][validation_size]:
					all_mets[model][validation_size][scenario] = cur_metrics[model]

	# print("\n=================================Final Result.")
	# print(json.dumps(all_mets))

	return all_mets


def plot_metric(all_mets, range_name, plot_dir):
	xs = range(1, 8)

	validation_size = [200, 400, 600, 800, 1000, 2000, 3000]

	# Create 5x3 subplots
	fig, axs = plt.subplots(5, 3, figsize=(10, 12))

	# Draw lines.
	
	ax_x = 0
	ax_y = 0
	for metric in ["auc", "accuracy", "precision", "recall", "f1"]:
		for scenario in ["bothin", "onlyone", "neither"]:
			for model in ["pdm", "nw", "sw", "blast"]:
				cur_metrics = []
				for x in xs:
					met = all_mets[model][validation_size[x-1]][scenario][metric]
					cur_metrics.append(met)
				
				axs[ax_x][ax_y].plot(xs, cur_metrics, '-o', label=model)
				axs[ax_x][ax_y].set_xticks(xs)
				axs[ax_x][ax_y].set_xticklabels(validation_size)
			
			axs[ax_x][ax_y].set_title(f"{scenario}, {metric}")
			# axs[ax_x][ax_y].set_xlabel('Validation Set Size.')
			if ax_x == 0 and ax_y == 0:
				axs[ax_x][ax_y].legend()
			
			ax_y += 1
		
		ax_x += 1
		ax_y = 0

	fig.suptitle(f"Avg metrics of {range_name}.", fontsize=16)

	# Show the plot
	plt.tight_layout()  # Adjust subplots to fit into figure area.
	plt.savefig(plot_dir)


def obtain_cr_metric_values(base_dir):
	all_mets = {}

	for validation_size in [200, 400, 600, 800, 1000, 2000, 3000]:
		for scenario in ["bothin", "onlyone", "neither"]:
			for model in ["cl_1000000", "cl_1000001", "cl_1000002", "cl_1000003", "cl_1000004"]:
				print(f"\n========== {validation_size} ========== {scenario} =========== {model} ===============")
				
				rs_dir = f"{base_dir}/{model}/result/{validation_size}/pdms/{scenario}"
				print(f"Dealing with {rs_dir}")

				if model not in all_mets:
					all_mets[model] = {}
				if validation_size not in all_mets[model]:
					all_mets[model][validation_size] = {}
				if scenario not in all_mets[model][validation_size]:
					all_mets[model][validation_size][scenario] = evaluate.display_avg_cr(rs_dir)

	# print("\n=================================Final Result.")
	# print(json.dumps(all_mets))

	return all_mets


def plot_cr_metric(all_mets, metric, plot_dir):
	xs = range(1, 8)

	validation_size = [200, 400, 600, 800, 1000, 2000, 3000]

	# Create 5x3 subplots
	fig, axs = plt.subplots(5, 3, figsize=(10, 12))

	# Draw lines.
	
	ax_x = 0
	ax_y = 0
	for cl in ["1", "2", "3", "4", "weighted avg"]:
		for scenario in ["bothin", "onlyone", "neither"]:
			for model in ["cl_1000000", "cl_1000001", "cl_1000002", "cl_1000003", "cl_1000004"]:
				cur_metrics = []
				for x in xs:
					met = all_mets[model][validation_size[x-1]][scenario][cl][metric]
					cur_metrics.append(met)
				
				axs[ax_x][ax_y].plot(xs, cur_metrics, '-o', label=model)
				axs[ax_x][ax_y].set_xticks(xs)
				axs[ax_x][ax_y].set_xticklabels(validation_size)
			
			axs[ax_x][ax_y].set_title(f"Distance Class={cl}, {scenario}")
			# axs[ax_x][ax_y].set_xlabel('Validation Set Size.')
			if ax_x == 0 and ax_y == 0:
				axs[ax_x][ax_y].legend()
			
			ax_y += 1
		
		ax_x += 1
		ax_y = 0

	fig.suptitle(f"The trend of avg value of {metric} along validation size.", fontsize=16)

	# Show the plot
	plt.tight_layout()  # Adjust subplots to fit into figure area.
	plt.savefig(plot_dir)


def plot_cr_accuracy(all_mets, plot_dir):
	xs = range(1, 8)
	metric = "accuracy"

	validation_size = [200, 400, 600, 800, 1000, 2000, 3000]

	# Create 1x3 subplots
	fig, axs = plt.subplots(1, 3, figsize=(10, 3))

	# Draw lines.
	
	ax_y = 0
	for scenario in ["bothin", "onlyone", "neither"]:
		for model in ["cl_1000000", "cl_1000001", "cl_1000002", "cl_1000003", "cl_1000004"]:
			cur_metrics = []
			for x in xs:
				met = all_mets[model][validation_size[x-1]][scenario][metric]
				cur_metrics.append(met)
			
			axs[ax_y].plot(xs, cur_metrics, '-o', label=model)
			axs[ax_y].set_xticks(xs)
			axs[ax_y].set_xticklabels(validation_size)
		
		axs[ax_y].set_title(f"{scenario}")
		if ax_y == 0:
			axs[ax_y].legend()
		
		ax_y += 1

	fig.suptitle(f"The trend of avg value of {metric} along validation size.", fontsize=16)

	# Show the plot
	plt.tight_layout()  # Adjust subplots to fit into figure area.
	plt.savefig(plot_dir)


if __name__ == "__main__":
	base_dir = "/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/evaluation_result"
	
	# subset_name = "cl_1000001"
	# plot_dir = "/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/pics"
	# all_mets = obtain_metric_values(base_dir, subset_name)
	# plot_metric(all_mets, "CL=1000001, 0-2 vs. 3-4", f"{plot_dir}/avg_metric_cl_1000001_2.png")

	plot_dir = "/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/pics"
	all_mets = obtain_cr_metric_values(base_dir)
	plot_cr_metric(all_mets, "precision", f"{plot_dir}/avg_5_cl_precision.png")
	plot_cr_metric(all_mets, "recall", f"{plot_dir}/avg_5_cl_recall.png")
	plot_cr_metric(all_mets, "f1-score", f"{plot_dir}/avg_5_cl_f1.png")
	plot_cr_accuracy(all_mets, f"{plot_dir}/avg_5_cl_accuracy.png")
