import numpy as np
from scipy.optimize import brute
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define sensor models
def ir3_func(x, a, b, c):
    return a / x + b * x + c

def ir4_func(x, x1, x2, off, l1, i1):
    conditions = [x <= x1, x1 < x, x >= x2]
    functions = [lambda x: x * l1 + off, x1 * l1 + off, lambda x: i1 / x + x1 * l1 + off - i1 / x2]
    return np.piecewise(x, conditions, functions)

def sonar1_func(x, a, b):
    return a * x + b

def sonar2_func(x, x1, off, l1):
    conditions = [x <= x1, x > x1]
    functions = [off, lambda x: (x - x1) * l1 + off]
    return np.piecewise(x, conditions, functions)

def inverse_sonar1_func(y, a, b):
    return (y - b) / a

# Inverse sensor model
def h_inverse(h, params, z, xmin, xmax):
    def f(x, z):
        return abs(h(x, *params) - z) ** 2

    xest = brute(f, ((xmin, xmax),), (z,))[0]
    return xest

# Motion model
def motion_model(position, velocity, acceleration, dt):
    return position + velocity * dt + (acceleration * dt**2)/2

# Calculate the process nosie variance
def process_noise_var(vel, dist, time):
    Wn = [dist[0]]

    for i in range(1, len(time)):
        dt = time[i] - time[i-1]
        acc = (vel[i]-vel[i-1])/dt
        Wn.append(dist[i] - motion_model(dist[i-1], vel[i-1], acc, dt))

    return np.var(Wn)

# Calculate the local slope of a function at a given point
def calculate_local_slope(func, params, x_best, epsilon=1e-6):
    slope = (func(x_best + epsilon, *params) - func(x_best, *params)) / epsilon

    if abs(slope) < 1e-1:
        return 1e-1
    else:
        return slope

# Get the variance of the estimates
def calculate_estimator_var(func, params, est, sensor_noise_var):
    local_slope = calculate_local_slope(func, params, est)
    sensor_estimate_var = 1/(local_slope**2) * sensor_noise_var

    if (sensor_estimate_var > 1000):
        return 1000
    else:
        return sensor_estimate_var

# Calculate the residual variance in mutiple windows
def calculate_residual_variances(x, y, model_func, params, length, num_windows=50):
    window_size = len(x) // num_windows
    mask_size = length // num_windows
    last_window_size = length - mask_size * num_windows
    residual_variances = []

    for i in range(num_windows):
        window_x = x[i * window_size: (i + 1) * window_size]
        window_y = y[i * window_size: (i + 1) * window_size]
        window_residuals = window_y - model_func(window_x, *params)
        window_residual_variance = np.var(window_residuals)
        residual_variances.extend([window_residual_variance] * mask_size)

    # Handle the last window
    if len(x) % num_windows != 0:
        last_window_x = x[num_windows * window_size:]
        last_window_y = y[num_windows * window_size:]
        last_window_residuals = last_window_y - model_func(last_window_x, *params)
        last_window_residual_variance = np.var(last_window_residuals)
        residual_variances.extend([last_window_residual_variance] * last_window_size)

    return residual_variances

# Function to find sensor model parameters and variances
def find_sensor_model_and_variances(x, y, model_func, length, threshold=0.025):
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

    residual_variances = calculate_residual_variances(x, y, model_func, params, length)

    return params, residual_variances

# Function to get the combined BLUE estimate and corresponding variance
def calculate_combined_estimate(*args):
    combined_est_num = 0
    combined_est_denom = 0
    combined_est_var_denom = 0

    for sensor_data, sensor_var in args:
        combined_est_num += sensor_data / sensor_var
        combined_est_denom += 1 / sensor_var
        combined_est_var_denom += 1 / sensor_var

    combined_est = combined_est_num / combined_est_denom
    combined_est_var = 1 / combined_est_var_denom

    return combined_est, combined_est_var

##### DATASET PREFERENCES #####
test_dataset = True # Change to true to perform EKF on the test dataset (must use training1 dataset for test)
training_dataset = 'training1.csv' # Either 'training1.csv' or 'training2.csv'
outlier_factor = 1000 # The degree to which outliers should be squashed

# Load training data
data = np.loadtxt(training_dataset, delimiter=',', skiprows=1)
index, time, distance, velocity_command, raw_ir1, raw_ir2, raw_ir3, raw_ir4, sonar1, sonar2 = data.T

# Get length of test dataset
if (test_dataset):
    data = np.loadtxt('test.csv', delimiter=',', skiprows=1)
    index, _, _, _, _, _, _, _, _ = data.T

# Find the models for each sensor
ir3_params, ir3_var = find_sensor_model_and_variances(distance, raw_ir3, ir3_func, len(index))
ir4_params, ir4_var = find_sensor_model_and_variances(distance, raw_ir4, ir4_func, len(index))
sonar1_params, sonar1_var = find_sensor_model_and_variances(distance, sonar1, sonar1_func, len(index))
# sonar2_params, sonar2_var = find_sensor_model_and_variances(distance, sonar2, sonar2_func, len(index))

# Load the rest of the test data
if (test_dataset):
    data = np.loadtxt('test.csv', delimiter=',', skiprows=1)
    index, time, velocity_command, raw_ir1, raw_ir2, raw_ir3, raw_ir4, sonar1, sonar2 = data.T

# Initial conditions
ir3_initial = h_inverse(ir3_func, ir3_params, raw_ir3[0], 0, 0.1)
ir4_initial = h_inverse(ir4_func, ir4_params, raw_ir4[0], 0, 0.1)
sonar1_initial = h_inverse(sonar1_func, sonar1_params, sonar1[0], 0, 0.1)
prev_position_prediction = (ir3_initial+ir4_initial+sonar1_initial)/3 # Average of initial sensor values
process_noise_var = process_noise_var(velocity_command, distance, time)
position_prediction_var = process_noise_var

# Track estimates and variances
position_estimates = [prev_position_prediction]
position_prediction_variances = [0]
sensor_weighting = [0]
model_weighting = [0]

# Loop through time steps
for i in range(1, len(index)):
    # Time difference between current and previous time step
    dt = time[i] - time[i-1]

    # Prediction Step
    acceleration = (velocity_command[i]-velocity_command[i-1])/dt
    position_prediction = motion_model(prev_position_prediction, velocity_command[i-1], acceleration, dt)
    position_prediction_var += process_noise_var # Update variance

    # Limits for h_inverse()
    xmin = position_prediction-0.1
    xmax = position_prediction+0.1

    # Sensor Measurement Step
    ir3_estimate = h_inverse(ir3_func, ir3_params, raw_ir3[i], xmin, xmax)
    ir3_estimate_var = calculate_estimator_var(ir3_func, ir3_params, prev_position_prediction, ir3_var[i])

    ir4_estimate = h_inverse(ir4_func, ir4_params, raw_ir4[i], xmin, xmax)
    ir4_estimate_var = calculate_estimator_var(ir4_func, ir4_params, prev_position_prediction, ir4_var[i])

    sonar1_estimate = inverse_sonar1_func(sonar1[i], *sonar1_params)
    sonar1_estimate_var = calculate_estimator_var(sonar1_func, sonar1_params, prev_position_prediction, sonar1_var[i])

    # sonar2_estimate = h_inverse(sonar2_func, sonar2_params, sonar2[i], xmin, xmax)
    # sonar2_estimate_var = calculate_estimator_var(position_estimates, dt, sonar2_var[i])

    # Ignore obvious outliers
    if ( abs(ir3_estimate-prev_position_prediction) > 0.1 ):
        ir3_estimate_var *= outlier_factor

    if ( abs(ir4_estimate-prev_position_prediction) > 0.1 ):
        ir4_estimate_var *= outlier_factor

    if ( abs(sonar1_estimate-prev_position_prediction) > 0.1 ):
        sonar1_estimate_var *= outlier_factor

    # Limit sensors within rated ranges
    argument_tuples = []
    if prev_position_prediction >= 0.1 and prev_position_prediction <= 0.8:
        argument_tuples.append((ir3_estimate, ir3_estimate_var))

    if prev_position_prediction >= 1 and prev_position_prediction <= 5:
        argument_tuples.append((ir4_estimate, ir4_estimate_var))

    if prev_position_prediction >= 0.02 and prev_position_prediction <= 4:
        argument_tuples.append((sonar1_estimate, sonar1_estimate_var))

    # if prev_position_prediction >= 0.3 and prev_position_prediction <= 5:
    #     argument_tuples.append((sonar2_estimate, sonar2_estimate_var))

    # Sensor Fusion (BLUE Estimator) 
    if (argument_tuples):
        combined_est, combined_est_var = calculate_combined_estimate(*argument_tuples)
    else:
        # Ignore sensors if the current best estimate is not within rated ranges
        combined_est, combined_est_var = position_prediction, position_prediction_var
        
    if ( abs(combined_est-prev_position_prediction) > 0.1 ):
        combined_est_var *= outlier_factor

    # Kalman Gain Calculation
    kalman_gain = 1/combined_est_var / (1/position_prediction_var + 1/combined_est_var)
    
    # Updated Position Prediction using Kalman Gain
    position_prediction = kalman_gain * combined_est + (1 - kalman_gain) * position_prediction

    # Update Variance of BLUE Estimator
    position_prediction_var = (combined_est_var*position_prediction_var)/(combined_est_var+position_prediction_var)

    # Store Predicted Position and Continue to Next Step
    prev_position_prediction = position_prediction
    position_estimates.append(position_prediction)
    position_prediction_variances.append(position_prediction_var)
    sensor_weighting.append(kalman_gain)
    model_weighting.append(1-kalman_gain)

# Position plot
pos_fig, pos_axes = plt.subplots(1)
if (not test_dataset):
    pos_axes.plot(time, distance, label='Actual position')
pos_axes.plot(time, position_estimates, label='Estimated position')
pos_axes.set_xlabel('Time (s)')
pos_axes.set_ylabel('Distance (m)')
pos_axes.legend()

# Weights plot
wgt_fig, wgt_axes = plt.subplots(1)
wgt_axes.plot(time, sensor_weighting, label='Estimate Weight')
wgt_axes.plot(time, model_weighting, label='Model Weight')
wgt_axes.set_xlabel('Time (s)')
wgt_axes.set_ylabel('Weight')
wgt_axes.legend()

# Show the plots
plt.show()