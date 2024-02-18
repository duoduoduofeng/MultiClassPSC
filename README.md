## Introduction
This project is similar to ProteinStructureClassification, but treat the problem as a multi classification.

## Architecture
- 📁 **data/** - Datasets downloaded from official SCOP website.
- 📁 **src/**
  - 📄 **main.py** - The entrance of PDM. Includes training and predicting.
  - 📄 **predict.py** - Use this for prediction. Avoid OOM with data size no greater than 800.
  - 📄 **batch_predict.py** - For batch predicting evaluation files.
  - 📄 **visualize_model.py** - Visualizes the PDM model structure.
  - 📁 **model/** - Contains the model architecture.
  - 📁 **preprocess/** - Constructs dataset for PDM.
  - 📁 **prepare_whole_dataset/** - How to build dataset from SCOP datasets.
  - 📁 **evaluates/** - Evaluates PDM and benchmarks.
- 📁 **generated_data/** - Stores generated intermediate data.
