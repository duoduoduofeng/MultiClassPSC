## Introduction
This project is similar to ProteinStructureClassification, but treat the problem as a multi classification.

## Architecture
- data
The datasets downloaded from official SCOP websit.
- src
	- main.py
	The start. It includes training and predicting.
	- predict.py
	If you only want to do prediction, use this script. Remember the target data size can not be larger than 800, otherwise OOM will occur.
	- batch_predict.py
	This is a script for batch predicting for the evaluation files.
	- visualize_model.py
	This script is used to visualize the inner structure of our PDM model.
	- model
	The model architecture.
	- preprocess
	Construct the designated dataset for PDM.
	- prepare_whole_dataset
	How do we build dataset from SCOP datasets.
	- evaluates
	How do we evaluate PDM and the benchmarks.
- generated_data
The generated intermediate data is stored here.