# Matlab PMML Demo
Simple application of MATLAB-PMML package for the prediction of energy consumption of a milling machine. 

## About
* Trains a Gaussian Process Regression Model to predict the energy consumption of a milling machine
* Saves the trained models to PMML using the (MATLAB PMML Package)[https://github.com/maxkferg/matlab-pmml]
* Makes energy predictions using the PMML model files

## Usage
Clone the repository and setup the dependancies
```sh
git clone https://github.com/maxkferg/matlab-pmml-demo
git submodule init
git submodule update
```

From the MATLAB editor:

1. Open the root folder with MATLAB
2. Run the train_energy_model.m file from MATLAB to generate PMML models
3. Run the predict_energy_usage.m file from MATLAB to make predictions

## License
MIT
