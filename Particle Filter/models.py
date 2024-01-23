"""Particle filter sensor and motion model implementations."""

import numpy as np
from numpy import cos, sin, tan, arccos, arcsin, arctan2, sqrt, exp
from numpy.random import randn
from utils import gauss, wraptopi, angle_difference


def motion_model(particle_poses, speed_command, odom_pose, odom_pose_prev, dt):
    """Apply motion model and return updated array of particle_poses.

    Parameters
    ----------

    particle_poses: an M x 3 array of particle_poses where M is the
    number of particles.  Each pose is (x, y, theta) where x and y are
    in metres and theta is in radians.

    speed_command: a two element array of the current commanded speed
    vector, (v, omega), where v is the forward speed in m/s and omega
    is the angular speed in rad/s.

    odom_pose: the current local odometry pose (x, y, theta).

    odom_pose_prev: the previous local odometry pose (x, y, theta).

    dt is the time step (s).

    Returns
    -------
    An M x 3 array of updated particle_poses.

    """

    M = particle_poses.shape[0]
    
    v = speed_command[0];
    w = speed_command[1];
    
    local_phi1 = arctan2((odom_pose[1]-odom_pose_prev[1]),(odom_pose[0]-odom_pose_prev[0]))-odom_pose_prev[2]
    local_d = ((odom_pose[1]-odom_pose_prev[1])**2 + (odom_pose[0]-odom_pose_prev[0])**2)**(1/2)
    local_phi2 = odom_pose[2] - odom_pose_prev[2] - local_phi1
    
    for m in range(M):
        particle_poses[m,0] = particle_poses[m,0] + local_d*cos(particle_poses[m,2]+local_phi1)
        particle_poses[m,1] = particle_poses[m,1] + local_d*sin(particle_poses[m,2]+local_phi1)
        particle_poses[m,2] = particle_poses[m,2] + local_phi1 + local_phi2
        
    for m in range(M):
        check = randn(1)
        value = 1
        if check < 0.5:
            value = -1
        particle_poses[m, 0] += randn(1)*0.1 * 0.4 *value#*randn(1)
        particle_poses[m, 1] += randn(1)*0.1 * 0.4 *value
        particle_poses[m, 2] += randn(1)*0.1 * 0.05 *value
    
    return particle_poses


def sensor_model(particle_poses, beacon_pose, beacon_loc):
    """Apply sensor model and return particle weights.

    Parameters
    ----------
    
    particle_poses: an M x 3 array of particle_poses (in the map
    coordinate system) where M is the number of particles.  Each pose
    is (x, y, theta) where x and y are in metres and theta is in
    radians.

    beacon_pose: the measured pose of the beacon (x, y, theta) in the
    robot's camera coordinate system. 

    beacon_loc: the pose of the currently visible beacon (x, y, theta)
    in the map coordinate system.

    Returns
    -------
    An M element array of particle weights.  The weights do not need to be
    normalised.

    """
 
    M = particle_poses.shape[0]
    particle_weights = np.zeros(M)
    
    r = sqrt((beacon_pose[0])**2 + (beacon_pose[1])**2)
    theta = arctan2(beacon_pose[1],beacon_pose[0]) 
    # For each particle calculate its weight based on its pose,
    # the relative beacon pose, and the beacon location.

    for m in range(M):
        r_m = sqrt((beacon_loc[0]-particle_poses[m,0])**2 + (beacon_loc[1]-particle_poses[m,1])**2)
        theta_m = angle_difference(particle_poses[m,2],arctan2(beacon_loc[1]-particle_poses[m,1],beacon_loc[0]-particle_poses[m,0]))
        
        particle_weights[m] = gauss(r-r_m,0,.1)*gauss(angle_difference(theta_m,theta),0,.075)
        #particle_weights[m] = 1;
    

    return particle_weights