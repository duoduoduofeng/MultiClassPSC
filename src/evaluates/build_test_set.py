import sys
import os
import random


def build_neither_validation_set(input_file, output_dir, \
	sample_size, sample_file_num, mode = "neither", theseed = 2024):
	file_row_count = 0
	
	# Count lines.
	with open(input_file, 'r') as fin:
		for line in fin:
			file_row_count += 1
	print(f"There are {file_row_count} rows in file {input_file}.\n")

	# Sampling one file by one file.
	# Here multipling 1.2 because the sampled size may be less than the target.
	prob = float(sample_size * 1.2 / file_row_count)
	for i in range(sample_file_num):
		output_file = f"{output_dir}/{mode}/{i+1}.txt"
		sample_volume = sample_size
		random.seed(theseed + i)
		with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
			row_count = 0
			sampled_row_count = 0
			for line in fin:
				row_count += 1
				if row_count == 1:
					fout.write(line)
				if sample_volume == 0:
					continue
				if random.random() < prob:
					fout.write(line)
					sampled_row_count += 1
					sample_volume -= 1

		print(f"Finished sampling {output_file} with {sampled_row_count} rows.")

	print(f"Finished sampling all output files.\n")


def load_trained_prots(train_set_file):
	row_count = 0
	prot_dict = {}
	with open(train_set_file, 'r') as fin:
		for line in fin:
			row_count += 1
			if row_count == 1:
				continue
			parts = line.strip().split("\t")
			if parts[0] not in prot_dict:
				prot_dict[parts[0]] = 1;
			if parts[1] not in prot_dict:
				prot_dict[parts[1]] = 1

	return prot_dict


def build_validation_set(input_file, output_dir, train_set_file, \
	sample_size, sample_file_num, mode = "bothin", theseed = 2024):

	# Load proteins in the training set.
	prot_dict = load_trained_prots(train_set_file)
	
	file_row_count = 0
	
	# Count lines.
	with open(input_file, 'r') as fin:
		for line in fin:
			file_row_count += 1
	print(f"There are {file_row_count} rows in file {input_file}.\n")

	# Sampling one file by one file.
	prob = float(sample_size * 4 / file_row_count)
	for i in range(sample_file_num):
		output_file = f"{output_dir}/{mode}/{i+1}.txt"
		sample_volume = sample_size
		random.seed(theseed + i)
		with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
			row_count = 0
			sampled_row_count = 0
			for line in fin:
				row_count += 1
				if row_count == 1:
					fout.write(line)
				if sample_volume == 0:
					continue
				if random.random() < prob:
					parts = line.strip().split("\t")
					if mode == "bothin":
						if parts[0] not in prot_dict or parts[1] not in prot_dict:
							continue
					elif mode == "onlyone":
						if parts[0] in prot_dict and parts[1] in prot_dict:
							continue
						elif parts[0] not in prot_dict and parts[1] not in prot_dict:
							continue
					
					fout.write(line)
					sampled_row_count += 1
					sample_volume -= 1

		print(f"Finished sampling {output_file} with {sampled_row_count} rows.")

	print(f"Finished sampling all output files.\n")



if __name__ == "__main__":
	base_dir = "/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data"
	model_name = "cl_1000004"
	
	train_set_file = f"{base_dir}/datasets/try_{model_name}/data/sample_proteins_dataset.train.txt"
	validation_file = f"{base_dir}/datasets/try_{model_name}/data/sample_proteins_dataset.validate.txt"
	excluded_validation_file = f"{base_dir}/datasets/try_{model_name}/data/sample_proteins_dataset.excluded_validate.txt"

	# sample_size = 200
	sample_file_num = 10
	for sample_size in [200, 400, 600, 800, 1000, 2000, 3000]:
		output_dir = f"{base_dir}/evaluation_result/{model_name}/testset/{sample_size}"
		build_validation_set(validation_file, output_dir, train_set_file, \
			sample_size, sample_file_num, mode = "bothin")
		build_validation_set(validation_file, output_dir, train_set_file, \
			sample_size, sample_file_num, mode = "onlyone")

		build_neither_validation_set(excluded_validation_file, output_dir, \
			sample_size, sample_file_num)


