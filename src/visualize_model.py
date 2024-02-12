import torch
import torch.nn as nn
import torch.optim as optim
from torch.nn.utils.rnn import pad_sequence
from torch.utils.data import DataLoader, random_split
from torchviz import make_dot

from model.protein_distance_model import *


def generate_toy_example():
    sequences1 = [
        "TVDVGPDSVKSACIEVDIQQTFFDKTWPRPIDVSKADGIIYPQGRT"
    ]
    sequences2 = [
        "QVYNFKRLVFTNCNYN"
    ]

    distances = [2]

    # Convert sequences to numerical indices
    abbrs = "ACDEFGHIKLMNPQRSTVWYX"
    char_to_index = {char: i for i, char in enumerate(abbrs)}
    sequences1 = [torch.tensor([char_to_index[char] for char in seq]) \
                        for seq in sequences1]
    sequences2 = [torch.tensor([char_to_index[char] for char in seq]) \
                        for seq in sequences2]
    # Pad sequences, which is necessary
    sequences1 = pad_sequence(sequences1, batch_first=True)
    sequences2 = pad_sequence(sequences2, batch_first=True)

    real_dis_tensor = torch.tensor(distances, dtype=torch.long).view(-1, 1)

    return sequences1, sequences2, real_dis_tensor


def vis_model(archi_file, the_device = "cpu", class_num = 5):
    the_embedding_dim = 128
    the_hidden_dim = 128
    
    # Set the model to evaluation mode (important if you have dropout layers)
    model = ProteinDistanceModel(
        embedding_dim=the_embedding_dim, 
        hidden_dim=the_hidden_dim, 
        num_classes=class_num)
    
    # Dom't load trained model.

    sequences1, sequences2, real_dis_tensor = generate_toy_example()
    predictions = model(sequences1, sequences2, real_dis_tensor)

    # print(dict(model.named_parameters()))
    print(model)
    
    # show_attrs=True, show_saved=True
    make_dot(predictions, 
        params=dict(list(model.named_parameters()))).\
        render(archi_file, format="png", cleanup=True)


if __name__ == "__main__":
    archi_file = "/Users/duoduo/Documents/lifeInCA/studyInTRU/2023Fall/graduate_project/other_psc/MultiClassPSC/generated_data/pics/PDM_classification"
    vis_model(archi_file)




