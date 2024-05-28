import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import Pipeline
from sklearn.metrics import mean_squared_error, r2_score
from nyoka import skl_to_pmml

# Set font properties for plots
rcParams['font.family'] = 'Times New Roman'
rcParams['font.size'] = 10

# Path to the CSV file using relative path
data_filename = "Drying fan speed.csv"  # Change to your desired file
data_path = os.path.join("data", data_filename)  # Assuming the data folder is at the root of the project
print(f"Data path: {data_path}")

# Load data
df = pd.read_csv(data_path, sep=';', decimal=",")
X = df[['Drying fan speed in rpm']]  # Change this to match the column in your CSV
y = df['P_avg in W']  # Change this if needed

def fit_score_plot_regression(model, X, y, xlabel, ylabel):
    """Fit model, plot regression, and display RMSE and R^2."""
    model.fit(X, y)
    y_predicted = model.predict(X)
    rmse = np.sqrt(mean_squared_error(y, y_predicted))
    r_value = r2_score(y, y_predicted)
    
    print(f"RMSE = {rmse:.2f}, R^2 = {r_value:.2f}")
    
    plt.scatter(X.values.ravel(), y, color="black", alpha=0.5, label="Sample data")
    plt.plot(X.values.ravel(), y_predicted, color="black", label="Regression model")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.legend()
    plt.show()

# Define the polynomial regression model with a pipeline
polynomial_regression = Pipeline([
    ('poly_features', PolynomialFeatures(degree=3, include_bias=False)),
    ('linear_regression', LinearRegression())
])

# Plot polynomial regression
fit_score_plot_regression(polynomial_regression, X, y, r'$\mathit{p}_{\mathrm{cleaning}}$ in bar', r'$\overline{P}$ in W' )

# Fit the model
polynomial_regression.fit(X, y)

# Export the trained model as PMML
pmml_file_path = os.path.join('models', 'model_cleaning_pump.pmml')  # Assuming 'models' directory at project root
skl_to_pmml(polynomial_regression, X.columns.tolist(), 'P_avg in W', pmml_file_path)
print(f'Model exported as PMML to: {pmml_file_path}')
