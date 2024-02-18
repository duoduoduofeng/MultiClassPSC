import os
import sys

# This script is used to split the mix validation file into neither, onlyone and both.

def get_uniq_pros(train_set_file):
	uniq_pros = {}

	with open(train_set_file, 'r') as fin:
		row_num = 0
		for line in fin:
			row_num += 1

			if row_num == 1:
				continue

			parts = line.strip().split("\t")

			if parts[0] not in uniq_pros:
				uniq_pros[parts[0]] = 0
			uniq_pros[parts[0]] += 1
			if parts[1] not in uniq_pros:
				uniq_pros[parts[1]] = 0
			uniq_pros[parts[1]] += 1

	print(f"There are {len(uniq_pros)} uniq proteins inside the training set.")
	return uniq_pros


def divide(uniq_pros, whole_set_file):
	bothin_file = f"{whole_set_file}.bothin"
	onlyone_file = f"{whole_set_file}.onlyone"
	neither_file = f"{whole_set_file}.neither"

	with open(whole_set_file, 'r') as fin, \
		open(bothin_file, 'w') as fout1, \
		open(onlyone_file, 'w') as fout2, \
		open(neither_file, 'w') as fout3:
		row_num = 0
		bothin_count = 0
		onlyone_count = 0
		neither_count = 0
		first_line = ""


		for line in fin:
			row_num += 1

			if row_num == 1:
				first_line = line
				fout1.write(first_line)
				fout2.write(first_line)
				fout3.write(first_line)
				continue

			parts = line.strip().split("\t")

			if parts[0] in uniq_pros and parts[1] in uniq_pros:
				bothin_count += 1
				fout1.write(line)
			elif parts[0] not in uniq_pros and parts[1] not in uniq_pros:
				neither_count += 1
				fout3.write(line)
			else:
				onlyone_count += 1
				fout2.write(line)

		print(f"Finished dividing file. bothin_count = {bothin_count}, onlyone_count = {onlyone_count}, neither_count = {neither_count}.")


if __name__ == "__main__":
	train_set_file = "/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/model_bak/MultiClassPSC/generated_data/datasets/try_cl_1000001.5class/data/sample_proteins_dataset.train.txt"
	whole_set_file = "/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/model_bak/MultiClassPSC/generated_data/datasets/try_cl_1000001.5class/data/sample_proteins_dataset.validate.txt.whole"
	
	uniq_pros = get_uniq_pros(train_set_file)
	divide(uniq_pros, whole_set_file)

