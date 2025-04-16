## Introduction
This project implements a Protein Distance Model (PDM) to predict distances between protein pairs according to the SCOP (Structural Classification of Proteins) standards. Combined with our aggregation strategy, it can ultimately predict the category of a single protein within SCOP. The primary objective of this model is to enhance the accuracy and precision of protein structure classification. Please contact newfdkl@gmail.com if you have any question.

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
- 📁 **results/** - Stores visual results of this project like charts, flowchart, etc.
