from datetime import datetime
import time

import torch
import torch.nn as nn
import torch.optim as optim
from torch.nn.utils.rnn import pad_sequence
from torch.utils.data import DataLoader, random_split

from model.process import *



if __name__ == "__main__":
    common_path = "../generated_data/datasets/try_mix"
    
    ### 1. Training
    dataset_file = f"{common_path}/data/sample_proteins_dataset.train.txt"
    epoch_times = 100

    timestamp = time.time()
    datetime_object = datetime.fromtimestamp(timestamp)
    dt_object = datetime_object.strftime('%Y%m%d_%H%M%S')

    the_device = "cuda:0" if torch.cuda.is_available() else "cpu"
    print("Using device:", the_device)

    model_save_file = f"{common_path}/models/trained_model.dict.{dt_object}"
    train_log = f"{common_path}/logs/train.log.{dt_object}"
    print(f"model_save_file: {model_save_file}, \ntrain_log: {train_log}\n")
    
    the_batch_size = 1500
    class_num = 6
    train(dataset_file, model_save_file, train_log, the_device, epoch_times, the_batch_size, class_num)

    ### 2. Prediction
    validate_dataset_file = f"{common_path}/data/sample_proteins_dataset.validate.txt"
    predict_result_file = f"{common_path}/result/predict_result.validate.txt.{dt_object}"
    predict(model_save_file, validate_dataset_file, predict_result_file, the_device, class_num)

    validate_dataset_file = f"{common_path}/data/sample_proteins_dataset.excluded_validate.txt"
    predict_result_file = f"{common_path}/result/predict_result.excluded_validate.txt.{dt_object}"
    predict(model_save_file, validate_dataset_file, predict_result_file, the_device, class_num)
