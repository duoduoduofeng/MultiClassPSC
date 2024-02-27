from datetime import datetime
import time

import torch
import torch.nn as nn
import torch.optim as optim
from torch.nn.utils.rnn import pad_sequence
from torch.utils.data import DataLoader, random_split

from model.process import *



if __name__ == "__main__":
    class_num = 5

    the_device = "cuda:0" if torch.cuda.is_available() else "cpu"
    print("Using device:", the_device)

    for i in range(10):
        file_no = i + 1
        
        trainset_dir = "try_cl_1000003"
        trained_model = "20240226_162531"
        
        validate_dataset_file = f"/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/evaluation_result/testset/neither/{file_no}.txt"
        predict_result_file = f"/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/evaluation_result/result/pdms/{trainset_dir}/{trained_model}/neither/{file_no}.txt"
        model_save_file = f"/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/datasets/{trainset_dir}/models/trained_model.dict.{trained_model}"
        predict(model_save_file, validate_dataset_file, predict_result_file, the_device, class_num)

        validate_dataset_file = f"/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/evaluation_result/testset/onlyone/{file_no}.txt"
        predict_result_file = f"/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/evaluation_result/result/pdms/{trainset_dir}/{trained_model}/onlyone/{file_no}.txt"
        predict(model_save_file, validate_dataset_file, predict_result_file, the_device, class_num)

        validate_dataset_file = f"/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/evaluation_result/testset/bothin/{file_no}.txt"
        predict_result_file = f"/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/evaluation_result/result/pdms/{trainset_dir}/{trained_model}/bothin/{file_no}.txt"
        predict(model_save_file, validate_dataset_file, predict_result_file, the_device, class_num)
