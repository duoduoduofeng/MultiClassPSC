from datetime import datetime
import torch
import torch.nn as nn
import torch.optim as optim
from torch.nn.utils.rnn import pad_sequence
from torch.utils.data import DataLoader, random_split

from preprocess.protein_dataset import *
from model.attention_layer import *
from model.protein_distance_model import *


def evaluation(model, data_loader, criterion):
    model.eval()
    total_loss = 0.0
    with torch.no_grad():
        for seq1, seq2, distance in data_loader:
            pred_distance = model(seq1, seq2, distance)
            loss = criterion(pred_distance, distance)
            total_loss += loss.item()

    avg_loss = total_loss / len(data_loader)
    return avg_loss


def train(dataset_file, model_save_file, train_log, the_device = 'cpu', \
        epoch_times = 10, the_batch_size = 256, class_num = 5):
    device = torch.device(the_device)
    # Some hyperparameters
    split_rate = 0.8

    learning_rate = 0.001

    # toy examples
    # the_batch_size = 2
    # epoch_times = 5
    
    # real dataset
    # the_batch_size = 128

    the_embedding_dim = 128
    the_hidden_dim = 128

    # toy example
    # sequences1, sequences2, distances = generate_toy_example()
    # And reflects the chars into numbers, by the operations in function load_dataset()

    # Load dataset
    print(f"Start loading dataset from {dataset_file}.\n")
    sequences1, sequences2, distances = load_dataset(dataset_file)
    print(f"Successfully loaded dataset with {len(sequences1)} sequences1, \
        {len(sequences2)} sequences2, \
        and {len(distances)} distances.\n")

    # Create dataset and DataLoader
    dataset = ProteinDataset(sequences1, sequences2, distances)
    train_size = int(split_rate * len(dataset))
    test_size = len(dataset) - train_size
    train_dataset, test_dataset = random_split(
        dataset, 
        [train_size, test_size])
    print(f"=************= Split datasets.\n")

    # Use pad_collate function for DataLoader
    train_loader = DataLoader(
        train_dataset, 
        batch_size=the_batch_size, 
        shuffle=True, 
        collate_fn=lambda batch: pad_collate(batch, device)
        )
    test_loader = DataLoader(
        test_dataset, 
        batch_size=the_batch_size, 
        shuffle=False, 
        collate_fn=lambda batch: pad_collate(batch, device)
        )
    print(f"=************= Load datasets, assigned pading function.\n")

    # Initialize the model
    model = ProteinDistanceModel(
        embedding_dim=the_embedding_dim, 
        hidden_dim=the_hidden_dim,
        num_classes = class_num).to(device)
    print(f"=************= Initialize model.\n")

    # Define loss function and optimizer
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    print(f"=************= Define loss function and optimizer.\n")

    # Training loop
    print(f"=************= Start training.........\n")
    with open(train_log, 'w') as fout:
        for epoch in range(epoch_times):
            model.train()
            for seq1, seq2, distance in train_loader:
                optimizer.zero_grad()
                pred_distance = model(seq1.to(device), seq2.to(device), distance.to(device))
                loss = criterion(pred_distance, distance)
                loss.backward()
                optimizer.step()

            # Evaluate on the test set
            avg_loss = evaluation(model, test_loader, criterion)
            ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            log_line = f"{ts}, Epoch {epoch + 1}, Test Loss: {avg_loss:.4f}"
            print(log_line)
            fout.write(f"{log_line}\n")


    print(f"\n=************= Finished training.\n")
    torch.save(model.state_dict(), model_save_file)
    print(f"Saved the trained model.\n")

    # return model



def predict(model_save_file, dataset_file, predict_result_file, the_device = "cpu", class_num = 5):
    print(f"=************= Start predicting... on {the_device}\n")

    the_embedding_dim = 128
    the_hidden_dim = 128
    
    # Set the model to evaluation mode (important if you have dropout layers)
    model = ProteinDistanceModel(
        embedding_dim=the_embedding_dim, 
        hidden_dim=the_hidden_dim, 
        num_classes=class_num)
    
    if the_device != "cpu":
        model.load_state_dict(torch.load(model_save_file))
        print(f"Loaded the trained model on gpu successfully.\n")
    else:
        model.load_state_dict(torch.load(model_save_file, map_location=torch.device(the_device)))
        print(f"Loaded the trained model on cpu successfully.\n")

    model.eval()

    sequences1, sequences2, real_dis_tensor, selected_test_set = load_test_data(dataset_file)
    # sequences1 = sequences1.to(the_device)
    # sequences2 = sequences2.to(the_device)
    # real_dis_tensor = real_dis_tensor.to(the_device)

    # Make predictions
    with torch.no_grad():
        predictions = model(sequences1, sequences2, real_dis_tensor)
        predict_distances = predictions.tolist()

    with open(predict_result_file, 'w') as fout:
        fout.write(f"protein1_pdb\tprotein2_pdb\tprotein1_classification\tprotein2_classfication\treal_distance\tpredict_distance\tpredict_class\n")
        for i in range(0, len(selected_test_set)):
            cur_record = selected_test_set[i]
            pro1 = cur_record["protein1_pdb"]
            pro2 = cur_record["protein2_pdb"]
            pro1_class = cur_record["protein1_classification"]
            pro2_class = cur_record["protein2_classfication"]
            real_distance = cur_record["distance"]
            predict_distance = predict_distances[i]
            predict_class = predict_distance.index(max(predict_distance))
            
            line = f"{pro1}\t{pro2}\t{pro1_class}\t{pro2_class}\t{real_distance}\t{predict_distance}\t{predict_class}"
            # print(line)
            fout.write(f"{line}\n")

    print(f"\n=************= Finished predicting.\n")
