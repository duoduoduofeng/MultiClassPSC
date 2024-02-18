from datetime import datetime
import time

import torch
import torch.nn as nn
import torch.optim as optim
from torch.nn.utils.rnn import pad_sequence
from torch.utils.data import DataLoader, random_split

from model.process import *



if __name__ == "__main__":
    dt_object = '20240131_154831'
    common_path = "../generated_data/datasets/try_cl_1000001"
    model_save_file = f"{common_path}/models/trained_model.dict.{dt_object}"
    class_num = 4

    the_device = "cuda:0" if torch.cuda.is_available() else "cpu"
    print("Using device:", the_device)

    validate_dataset_file = f"{common_path}/data/sample_proteins_dataset.validate.txt"
    predict_result_file = f"{common_path}/result/predict_result.validate.txt.{dt_object}"
    predict(model_save_file, validate_dataset_file, predict_result_file, the_device, class_num)

    validate_dataset_file = f"{common_path}/data/sample_proteins_dataset.excluded_validate.txt"
    predict_result_file = f"{common_path}/result/predict_result.excluded_validate.txt.{dt_object}"
    predict(model_save_file, validate_dataset_file, predict_result_file, the_device, class_num)
