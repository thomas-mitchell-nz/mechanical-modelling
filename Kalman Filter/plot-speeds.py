import numpy as np
import matplotlib.pyplot as plt

# Motion model func
def motion_model_func(position, velocity, acceleration, dt):
    return position + velocity * dt + (acceleration * dt**2)/2

# Returns motion model list and process nosie variance
def motion_model(vel, dist, time):
    prediction_list = [dist[0]]
    Wn = [0]

    for i in range(1, len(time)):
        dt = time[i] - time[i-1]
        acc = (vel[i] - vel[i-1]) / dt
        prediction = motion_model_func(dist[i-1], vel[i-1], acc, dt)
        prediction_list.append(prediction)
        Wn.append(dist[i] - prediction) 

    return prediction_list, np.var(Wn)

# Plot the speed data
def load_and_plot_data(filename, axis):
    # Load data
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    
    # Unpack data columns
    _, time, distance, v_com, _, _, _, _, _, _ = data.T

    # Get the position prediction and process noise variance from model
    position, wn_variance = motion_model(v_com, distance, time)
    print("Variance of Wn:", wn_variance)

    # Calculate speed predicted using your motion model
    v_model = np.gradient(position, time)

    # Plot speeds
    axis.plot(time, v_model, linewidth=1, label='Predicted speed')
    axis.plot(time, v_com, label='Commanded Speed')
    axis.set_xlabel('Time (s)')
    axis.set_ylabel('Speed (m/s)')
    axis.legend()
    
fig, axes = plt.subplots(1, 2)

# Load and plot data from training1.csv
load_and_plot_data('training1.csv', axes[0])

# Load and plot data from training2.csv
load_and_plot_data('training2.csv', axes[1])

plt.show()