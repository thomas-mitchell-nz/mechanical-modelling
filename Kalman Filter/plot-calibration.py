import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def positive_linear(x, a, b):
    return a * x + b

def inverse_model(x, a, b, c):
    return a / x + b * x + c

def piecewise_inverse(x, x1, x2, off, l1, i1):
    conditions = [x <= x1, x1 < x, x >= x2]
    functions = [lambda x: x * l1 + off, x1 * l1 + off, lambda x: i1 / x + x1 * l1 + off - i1 / x2]
    return np.piecewise(x, conditions, functions)

def piecewise_linear(x, x1, off, l1):
    conditions = [x <= x1, x > x1]
    functions = [off, lambda x: (x - x1) * l1 + off]
    return np.piecewise(x, conditions, functions)

def find_best_model(x, y, model_func, threshold=0.025):
    x_copy = np.array(x)
    y_copy = np.array(y)
    outliers_x = []
    outliers_y = []

    for _ in range(1000):
        params, _ = curve_fit(model_func, x_copy, y_copy)
        y_fit = model_func(x_copy, *params)
        residuals = y_copy - y_fit
        std_residuals = np.std(residuals)

        if std_residuals < threshold:
            break

        outlier_indices = np.abs(residuals) > 3 * std_residuals
        outliers_x.extend(x_copy[outlier_indices])
        outliers_y.extend(y_copy[outlier_indices])
        y_copy = y_copy[~outlier_indices]
        x_copy = x_copy[~outlier_indices]

    model_data = model_func(x_copy, *params)

    return x_copy, model_data, outliers_x, outliers_y, residuals

def plot_residuals(axis, x, y, label):
    axis.plot(x, y, '.', alpha=0.2, label=label + ' Residuals')
    axis.axhline(0, color='r', linestyle='--', linewidth=1)
    axis.set_xlabel('Distance (m)')
    axis.set_ylabel('Residuals')
    axis.set_title(label + ' Residuals')
    axis.legend()

def plot_data(axis, err_axis, x, y, label, model_func):
    x_model, model, x_outliers, y_outliers, residuals = find_best_model(x, y, model_func)
    axis.plot(x, y, '.', alpha=0.2, label=label)
    axis.plot(x_model, model, 'r-', label='Model')
    axis.plot(x_outliers, y_outliers, '.', label='Outliers')
    axis.set_xlabel('Distance (m)')
    axis.set_ylabel(label)
    axis.set_title(label)
    axis.legend()

    # Print model coefficients
    params, _ = curve_fit(model_func, x_model, model)
    print(f'{label} Coefficients:', params)

    # Calculate and print variance of residuals
    residuals_variance = np.var(residuals)
    print(f'{label} Residuals Variance:', residuals_variance)

    # Plot residuals on the error subplot
    plot_residuals(err_axis, x_model, residuals, label)

# Load data
filename = 'training1.csv'
data = np.loadtxt(filename, delimiter=',', skiprows=1)
index, time, distance, velocity_command, raw_ir1, raw_ir2, raw_ir3, raw_ir4, sonar1, sonar2 = data.T

# Model plot config
cal, cal_axes = plt.subplots(2, 3)
cal.suptitle('Calibration data')

# Error plot config
err, err_axes = plt.subplots(2, 3)
err.suptitle('Error data')

# Plot the data for each sensor
plot_data(cal_axes[0, 0], err_axes[0, 0], distance, raw_ir1, 'IR1', inverse_model)
plot_data(cal_axes[0, 1], err_axes[0, 1], distance, raw_ir2, 'IR2', inverse_model)
plot_data(cal_axes[0, 2], err_axes[0, 2], distance, raw_ir3, 'IR3', inverse_model)
plot_data(cal_axes[1, 0], err_axes[1, 0], distance, raw_ir4, 'IR4', piecewise_inverse)
plot_data(cal_axes[1, 1], err_axes[1, 1], distance, sonar1, 'Sonar1', piecewise_linear)
plot_data(cal_axes[1, 2], err_axes[1, 2], distance, sonar2, 'Sonar2', piecewise_linear)

# Show the plots
plt.show()