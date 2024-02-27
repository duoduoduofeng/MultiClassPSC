import sys
import json
import random


def extract_all_keys(train_file):
	all_trains = {}
	with open(train_file, 'r') as fin:
		row_num = 0
		for line in fin:
			row_num += 1
			if row_num < 2:
				continue

			parts = line.strip().split("\t")

			if parts[0] not in all_trains:
				all_trains[parts[0]] = {}
			if int(parts[2]) not in all_trains[parts[0]]:
				all_trains[parts[0]][int(parts[2])] = 0
			all_trains[parts[0]][int(parts[2])] += 0

			if parts[1] not in all_trains:
				all_trains[parts[1]] = {}
			if int(parts[2]) not in all_trains[parts[1]]:
				all_trains[parts[1]][int(parts[2])] = 0
			all_trains[parts[1]][int(parts[2])] += 0

	print(all_trains)
	return all_trains


def print_valid_pros(all_trains):
	valid_pros = {}
	for pro in all_trains:
		# for num in [0, 1, 2, 3, 4]:
		# 	if num not in all_trains[pro]:
		# 		break
		# 	valid_pros[pro] = all_trains[pro]
		thed = all_trains[pro]
		# if 0 in thed and 1 in thed and 2 in thed and 3 in thed and 4 in thed:
		if 0 in thed and 1 in thed and 2 in thed and 4 in thed and 8 in thed:
			valid_pros[pro] = all_trains[pro]

	print(valid_pros)
	return valid_pros


def show(train_file, valid_pros):
	has_print = {}
	to_print = {}

	with open(train_file, 'r') as fin:
		row_num = 0
		for line in fin:
			row_num += 1
			if row_num < 2:
				continue

			parts = line.strip().split("\t")

			if parts[0] in valid_pros:
				if parts[0] not in has_print:
					has_print[parts[0]] = []
				if int(parts[2]) not in has_print[parts[0]]:
					if parts[0] not in to_print:
						to_print[parts[0]] = []
					to_print[parts[0]].append(line.strip())
					has_print[parts[0]].append(int(parts[2]))
			elif parts[1] in valid_pros:
				if parts[1] not in has_print:
					has_print[parts[1]] = []
				if int(parts[2]) not in has_print[parts[1]]:
					if parts[1] not in to_print:
						to_print[parts[1]] = []
					to_print[parts[1]].append(line.strip())
					has_print[parts[1]].append(int(parts[2]))
	
	for key in to_print:
		print()
		print(key)
		for line in to_print[key]:
			print(line)



if __name__ == "__main__":
	train_file = sys.argv[1]
	all_trains = extract_all_keys(train_file)
	valid_pros = print_valid_pros(all_trains)
	show(train_file, valid_pros)
