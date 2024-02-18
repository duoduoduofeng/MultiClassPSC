## Introduction
This project is similar to ProteinStructureClassification, but treat the problem as a multi classification.

## Architecture
.
├── data
│   └── Datasets downloaded from official SCOP website.
├── src
│   ├── main.py
│   │   └── The start. It includes training and predicting.
│   ├── predict.py
│   │   └── Use this for prediction. Target data size cannot exceed 800 to avoid OOM.
│   ├── batch_predict.py
│   │   └── Script for batch predicting evaluation files.
│   ├── visualize_model.py
│   │   └── Visualizes the inner structure of the PDM model.
│   ├── model
│   │   └── Contains the model architecture.
│   ├── preprocess
│   │   └── Constructs the designated dataset for PDM.
│   ├── prepare_whole_dataset
│   │   └── How to build dataset from SCOP datasets.
│   └── evaluates
│       └── Evaluates PDM and benchmarks.
└── generated_data
    └── Stores generated intermediate data.
