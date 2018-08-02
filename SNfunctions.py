import numpy as np
from numpy import unravel_index
import math
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import copy

"""
Functions for locating the SN.
"""

# Function that takes in the latitude and longitude (in degrees) and returns the (x,y,z) coordinates
# of the point on earth (in m), with the centre of the earth as the origin
def find_xyz_from_latlong(latitude, longitude, radius):
    x = radius*math.cos(math.radians(longitude))*math.cos(math.radians(latitude))
    y = radius*math.sin(math.radians(longitude))*math.cos(math.radians(latitude))
    z = radius*math.sin(math.radians(latitude))
    return np.array([x,y,z])

# Function to display the latitude and longitude of points in 3D
def lat_long_from_xyz(point):
    lat = math.degrees(math.asin(point[2] / np.linalg.norm(point)))
    # If the point is on the north or south pole, the longitude can take any value
    # The following code makes sure the correct atan value is chosen, since in the range
    # -pi-->pi there are 2 atan values for any number
    if point[1] >= 0.:
        if point[0] >= 0.:
            long = math.degrees(math.atan(point[1]/point[0]))
        elif point[0] < 0.:
            long = 180. + math.degrees(math.atan(point[1]/point[0]))
    elif point[1] < 0.:
        if point[0] >= 0.:
            long = math.degrees(math.atan(point[1]/point[0]))
        elif point[0] < 0.:
            long = math.degrees(math.atan(point[1]/point[0])) - 180.
    return np.array([lat,long])

# Function that takes in 2 3D points and returns the separation
def distance_between_3d_points(p1, p2):
    separation_vector = p1 - p2
    return math.sqrt(separation_vector[0]**2+separation_vector[1]**2+separation_vector[2]**2)

# Function to produce the equation of the plane on which the SN lies
def produce_plane(pos1, pos2, dT, SN_dist, c_light):
    dT_max = distance_between_3d_points(pos1, pos2)/c_light
    theta = math.acos(dT/dT_max) # Find angle between the SN position and the relative positions of the detectors
    plane_normal = pos2-pos1
    mag_plane_normal = np.linalg.norm(plane_normal)
    plane_normal /= mag_plane_normal # Normalize the plane normal vector
    plane_point = SN_dist*math.cos(theta)*plane_normal
    plane = np.append(plane_normal, np.array([np.dot(plane_normal, plane_point)]))
    return plane

# Function to output the equation of a line intersecting 2 planes
def plane_intersection(plane_1, plane_2):
    # Define the line direction vector for planes 1 and 2
    line_vector = np.cross(plane_1[0:3], plane_2[0:3])
    # Initialise the point on the line
    line_point = np.empty(3)
    # Choose a valid point on the line
    if plane_1[0] != 0 and (plane_2[2]*plane_1[0]-plane_2[0]*plane_1[2]) != 0:
        # Define a point on the line (valid for plane_1[0] =/= 0)
        line_point[1] = 0.
        line_point[2] = (plane_1[0]*plane_2[3]-plane_2[0]*plane_1[3])/(plane_2[2]*plane_1[0]-plane_2[0]*plane_1[2])
        line_point[0] = (plane_1[3]/plane_1[0])-(plane_1[2]/plane_1[0])*line_point[2]
    elif plane_1[1] != 0 and (plane_2[1]*plane_1[2]-plane_2[2]*plane_1[1]) != 0:
        # Define a point on the line (valid for plane_1[1] =/= 0)
        line_point[0] = 0.
        line_point[1] = (plane_1[2]*plane_2[3]-plane_2[2]*plane_1[3])/(plane_2[1]*plane_1[2]-plane_2[2]*plane_1[1])
        line_point[2] = (plane_1[3]/plane_1[2])-(plane_1[1]/plane_1[2])*line_point[1]
    elif plane_1[2] != 0 and (plane_2[0]*plane_1[1]-plane_2[1]*plane_1[0]) != 0:
        # Define a point on the line (valid for plane_1[2] =/= 0)
        line_point[2] = 0.
        line_point[0] = (plane_1[1]*plane_2[3]-plane_2[1]*plane_1[3])/(plane_2[0]*plane_1[1]-plane_2[1]*plane_1[0])
        line_point[1] = (plane_1[3] / plane_1[1]) - (plane_1[0] / plane_1[1]) * line_point[0]
    # Return equation of line in form [gradx, pointx, grady, pointy, gradz, pointz]
    return np.array([line_vector[0],line_point[0],line_vector[1],line_point[1],line_vector[2],line_point[2]])

# Function to find the intersection points of 2 planes with a sphere
def intersection_2planes_sphere(line, sphere):
    a = line[0]**2 + line[2]**2 + line[4]**2
    b = 2*(line[0]*(line[1]-sphere[0])+line[2]*(line[3]-sphere[1])+line[4]*(line[5]-sphere[2]))
    c = (line[1]-sphere[0])**2 + (line[3]-sphere[1])**2 + (line[5]-sphere[2])**2 - sphere[3]**2
    # The check a==0 is to see if the planes are parallel, in which case they cannot intersect
    if a == 0:
        raise ValueError('Planes do not meet.')
    else:
        discriminant = b**2 - 4*a*c
        if discriminant >= 0:
            # Find the t value(s)
            t1 = (-b+math.sqrt(discriminant))/(2*a)
            t2 = (-b-math.sqrt(discriminant))/(2*a)
            # Find the positions of the intersection points by plugging the t value(s) into the line equation
            point1 = np.array([line[0]*t1+line[1],line[2]*t1+line[3],line[4]*t1+line[5]])
            point2 = np.array([line[0]*t2+line[1],line[2]*t2+line[3],line[4]*t2+line[5]])
            return point1, point2
        elif discriminant < 0:
            raise ValueError('Discriminant less than 0.')

def find_deltaT(pos1, pos2, posSN, c_light):
    local_posSN = copy.copy(posSN) # Define a local version of the SN position
    dT_max = distance_between_3d_points(pos1, pos2)/c_light
    relative_pos = pos2-pos1
    relative_pos /= np.linalg.norm(relative_pos)
    local_posSN /= np.linalg.norm(local_posSN)
    return (dT_max*np.dot(local_posSN, relative_pos))/(np.linalg.norm(local_posSN)*np.linalg.norm(relative_pos))

# Function that takes in quantities defining a circle in 3D space and returns a point on the circle at
# a user-given angle
# Theta is an angle in radians, centre and normal are 3-dimensional numpy arrays and radius is a number
def point_on_3d_circle(normal, centre, radius, theta):
    local_normal = copy.copy(normal) # Define a local variable for the normal vector so it is not changed outwith the function
    local_normal /= np.linalg.norm(local_normal) # normalize the normal vector
    # Now we find the first vector orthogonal to the normal vector
    # We do this by solving a.n = 0
    # If the normal vector is 0 then we cannot get a definite solution
    if local_normal[0] == local_normal[1] == local_normal[2] == 0.:
        raise ValueError('local_normal vector is the 0 vector.')
    if local_normal[0] == 0:
        a = np.array([1,0,0])
    elif local_normal[1] == 0.:
        a = np.array([0,1,0])
    elif local_normal[2] == 0.:
        a = np.array([0,0,1])
    else:
        a = np.array([1,1,(-local_normal[0]-local_normal[1])/local_normal[2]])
        a /= np.linalg.norm(a)
    # Our final orthogonal vector is found by b = n x a so it is perpendicular to both a and n
    b = np.cross(a,local_normal)
    # Now we define the x,y,z coordinates of the point on the circle using this method:
    # https://math.stackexchange.com/questions/73237/parametric-equation-of-a-circle-in-3d-space
    return centre+radius*math.cos(theta)*a+radius*math.sin(theta)*b

# Method to check if the SN is on the circle, taking the sphere radius and plane the circle lies on as input
def check_point_on_circle(plane, radius, point):
    lhs_plane = plane[0]*point[0]+plane[1]*point[1]+plane[2]*point[2]
    print "LHS plane equation = %f.  RHS plane equation = %f." %(lhs_plane, plane[3])
    print "Magnitude of point = %f.  Radius of sphere = %f" %(np.linalg.norm(point), radius)
    if lhs_plane == plane[3] and radius == np.linalg.norm(point):
        return True
    else:
        return False

# Function to produce a list of longitudes and latitudes of a 3D circle based on
# the ToA difference between 2 detectors
def circle_plot_points(N_points, pos1, pos2, c_light, dT, SN_dist):
    angle_list = np.linspace(0,2*np.pi,N_points)
    # Circle parameters
    theta = math.acos(dT*c_light/distance_between_3d_points(pos1, pos2))
    radius = SN_dist*math.sin(theta)
    circle_normal = pos2-pos1
    circle_normal /= np.linalg.norm(circle_normal) # Normalize the plane normal vector
    circle_centre = SN_dist*math.cos(theta)*circle_normal
    long_list = np.empty(N_points)
    lat_list = np.empty(N_points)
    # Fill arrays
    for i in range(len(angle_list)):
        angle = angle_list[i]
        current_point = point_on_3d_circle(circle_normal, circle_centre, radius, angle)
        current_lat_long = lat_long_from_xyz(current_point)
        lat_list[i] = current_lat_long[0]
        long_list[i] = current_lat_long[1]
    lat_list = lat_list
    long_list = long_list
    return lat_list, long_list

# Calculate the radius of a circle in the sky
def circle_radius(dT, c_light, pos1, pos2, SN_dist):
    theta = math.acos(dT * c_light / distance_between_3d_points(pos1, pos2))
    radius = SN_dist * math.sin(theta)
    return radius

# Gaussian distribution function adapted for the purposes of the problem
# Returns the floor of the function f(x) at a particular x
def gaussian_function(scale, mean_long, mean_lat, stddev, lat, long):
    f_lat = scale * math.exp((-(lat - mean_lat) ** 2)/(2 * stddev))
    f_long = scale * math.exp((-(long - mean_long) ** 2) / (2 * stddev))
    return math.floor((f_lat + f_long)/2.)

# Function to produce arrays of latitudes, longitudes and the corresponding weights for the 2D histogram
def produce_hist_arrays(lat_list, long_list, scale, circle_width, N_points, stddev):
    # Make sure the longitude and latitude lists are the same length
    if len(lat_list) != len(long_list):
        raise ValueError('The latitude list is not the same length as the longitude list.')
    # Initialise empty arrays to hold all latitudes, longitudes and their corresponding weights
    final_lat_list = np.empty([len(long_list), N_points])
    final_long_list = np.empty([len(long_list), N_points])
    weight_list = np.empty([len(long_list), N_points])
    # Loop through all points on the circle
    for i in range(len(lat_list)):
        # Retrieve coordinates
        current_lat = lat_list[i]
        current_long = long_list[i]
        # Create arrays of points surrounding the current point on the circle
        current_lat_list = np.empty(N_points)
        current_long_list = np.empty(N_points)
        # Find next lat/long in list
        if i != (len(lat_list)-1):
            next_lat = lat_list[i+1]
            next_long = long_list[i+1]
        else:
            next_lat = lat_list[0]
            next_long = long_list[0]
        # (lat,long) 2D coord of current circle point
        pos_current_point = np.array([current_lat, current_long])
        # Define normal vector to direction of travel of circle
        # If vector is (x,y) then (y,-x) is perpendicular to (x,y) in 2D coords
        normal_vector = np.array([next_long-current_long, current_lat-next_lat])
        # Normalize vector
        if np.linalg.norm(normal_vector) != 0.:
            normal_vector /= np.linalg.norm(normal_vector)
        scale_factor_list = np.linspace(-circle_width, circle_width, N_points)
        # Create array of points normal to the current direction of the circle at the location of the current point
        for k in range(N_points):
            temp_point = pos_current_point + scale_factor_list[k]*normal_vector
            current_lat_list[k] = temp_point[0]
            current_long_list[k] = temp_point[1]
        current_weight_list = np.empty(N_points)
        # Loop through array of surrounding points and get the weights
        for j in range(N_points):
            current_weight_list[j] = gaussian_function(scale, current_long, current_lat, stddev,
                                               current_lat_list[j], current_long_list[j])
        # Make sure the latitudes and longitude values calculated wrap round if they are outside the
        # allowed range - not sure it's entirely correct but it doesn't seem to cause problems
        for l in range(N_points):
            if current_lat_list[l] > 90.:
                current_lat_list[l] = 180. - current_lat_list[l]
            elif current_lat_list[l] < -90.:
                current_lat_list[l] = 180. + current_lat_list[l]
        for m in range(N_points):
            if current_long_list[m] > 180.:
                current_long_list[m] = current_long_list[m] - 360.
            elif current_long_list[m] < -180.:
                current_long_list[m] = current_long_list[m] + 360.
        # Add lat/longs to combined array
        final_lat_list[i] = current_lat_list
        final_long_list[i] = current_long_list
        weight_list[i] = current_weight_list
    # Reduce all lists to a 1D array
    final_long_list = final_long_list.flatten()
    final_lat_list = final_lat_list.flatten()
    weight_list = weight_list.flatten()
    return final_lat_list, final_long_list, weight_list

# Alex's method for finding the delta T values between detectors
def find_deltaT_alex(pos1, pos2, posSN, c_light):
    pos_1_SN = pos1 - posSN
    pos_2_SN = pos2 - posSN
    abs_1_SN = np.linalg.norm(pos_1_SN)
    abs_2_SN = np.linalg.norm(pos_2_SN)
    dr_1_2 = abs_1_SN-abs_2_SN
    dt_1_2 = dr_1_2/c_light
    return dt_1_2

# Function to produce a Skycoord object of long/lats of a 3D circle based on
# the ToA difference between 2 detectors
# The wrapped long/lats are returned
def circle_plot_points_astropy(N_points, pos1, pos2, c_light, dT, SN_dist):
    angle_list = np.linspace(0,2*np.pi,N_points)
    # Circle parameters
    theta = math.acos(dT*c_light/distance_between_3d_points(pos1, pos2))
    radius = SN_dist*math.sin(theta)
    circle_normal = pos2-pos1
    circle_normal /= np.linalg.norm(circle_normal) # Normalize the plane normal vector
    circle_centre = SN_dist*math.cos(theta)*circle_normal
    long_list = np.empty(N_points)
    lat_list = np.empty(N_points)
    # Fill arrays
    for i in range(len(angle_list)):
        angle = angle_list[i]
        current_point = point_on_3d_circle(circle_normal, circle_centre, radius, angle)
        current_lat_long = lat_long_from_xyz(current_point)
        lat_list[i] = current_lat_long[0]
        long_list[i] = current_lat_long[1]
    lat_list = lat_list
    long_list = long_list
    lat_list = lat_list * u.degree
    long_list = long_list * u.degree
    c = SkyCoord(ra=long_list, dec=lat_list, frame='icrs')
    ra_rad = c.ra.wrap_at(180 * u.deg).radian
    dec_rad = c.dec.radian
    return ra_rad, dec_rad
