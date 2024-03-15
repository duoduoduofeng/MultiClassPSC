from datetime import datetime
import time

import torch
import torch.nn as nn
import torch.optim as optim
from torch.nn.utils.rnn import pad_sequence
from torch.utils.data import DataLoader, random_split

from model.process import *



if __name__ == "__main__":
    the_device = "cuda:0" if torch.cuda.is_available() else "cpu"
    print("Using device:", the_device)

    base_dir = "/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data"
    
    trainset_dir = "cl_1000004"
    trained_model = "20240305_133008"
    model_save_file = f"{base_dir}/datasets/try_{trainset_dir}/models/trained_model.dict.{trained_model}"
    
    validation_size = 3000
    for scenario in ["bothin", "onlyone", "neither"]:
        for i in range(10):
            file_no = i + 1
            
            validate_dataset_file = f"{base_dir}/evaluation_result/{trainset_dir}/testset/{validation_size}/{scenario}/{file_no}.txt"
            predict_result_file = f"{base_dir}/evaluation_result/{trainset_dir}/result/{validation_size}/pdms/{scenario}/{file_no}.txt"
            print(f"=============== Dealing file {validate_dataset_file}.")
            predict(model_save_file, validate_dataset_file, predict_result_file, the_device, class_num = 5)