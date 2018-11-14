import copy
import math

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord


################################
# Functions for locating the SN.
################################

# Function that takes in the latitude and longitude (in degrees) and returns the (x,y,z) coordinates
# of the point on earth (in m), with the centre of the earth as the origin
#
# Input: latitude (number), longitude (number), radius (number)
#
# Output: numpy array [x, y, z]
#
def find_xyz_from_latlong(latitude, longitude, radius):
    x = radius*math.cos(math.radians(longitude))*math.cos(math.radians(latitude))
    y = radius*math.sin(math.radians(longitude))*math.cos(math.radians(latitude))
    z = radius*math.sin(math.radians(latitude))
    return np.array([x,y,z])

# Function that returns the latitude and longitude of a point in 3D
#
# Input: list/numpy array [x, y, z]
#
# Output: numpy array [latitude, longitude]
#
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
#
# Input: list/numpy array [x1, y1, z1], list/numpy array [x2, y2, z2]
#
# Output: number
#
def distance_between_3d_points(p1, p2):
    separation_vector = p1 - p2
    return math.sqrt(separation_vector[0]**2+separation_vector[1]**2+separation_vector[2]**2)

# Function to produce the equation of the plane on which the SN lies
#
# Input: list/numpy array [x1, y1, z1], list/numpy array [x2, y2, z2], number, number, number
#
# Output: numpy array [a, b, c, d]
#
def produce_plane(pos1, pos2, dT, SN_dist, c_light):
    dT_max = distance_between_3d_points(pos1, pos2) / c_light

    # Find angle between the SN position and the relative positions of the detectors
    theta = math.acos(dT / dT_max)

    plane_normal = pos2 - pos1
    mag_plane_normal = np.linalg.norm(plane_normal)
    plane_normal /= mag_plane_normal # Normalize the plane normal vector

    plane_point = SN_dist*math.cos(theta)*plane_normal

    # The plane equation is of the form [a, b, c, d] which defines a scalar equation of the plane
    # ax + by + cz = d
    plane = np.append(plane_normal, np.array([np.dot(plane_normal, plane_point)]))

    return plane

# Function to output the equation of a line intersecting 2 planes
#
# Input:  list/numpy array [a1, b1, c1, d1], list/numpy array [a2, b2, c2, d2]
#
# Output: numpy array [gradx, x, grady, y, gradz, z]
#
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

    # Return equation of line in form [gradx, x, grady, y, gradz, z]
    # [gradx, grady, gradz] refers to the gradient of the line
    # [x, y, z] refers to a point on the line
    return np.array([line_vector[0],line_point[0],line_vector[1],line_point[1],line_vector[2],line_point[2]])

# Function to find the intersection points of 2 planes with a sphere
# A sphere is defined as [x, y, z, r] where [x, y, z] is the centre and r is the radius
#
# Input: list/numpy array [gradx, x, grady, y, gradz, z], list/numpy array [x, y, z, r]
#
# Output: list/numpy array [x1, y1, z1], list/numpy array [x2, y2, z2]
#
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
            t1 = (-b+math.sqrt(discriminant)) / (2*a)
            t2 = (-b-math.sqrt(discriminant)) / (2*a)
            # Find the positions of the intersection points by plugging the t value(s) into the line equation
            point1 = np.array([line[0]*t1+line[1], line[2]*t1+line[3], line[4]*t1+line[5]])
            point2 = np.array([line[0]*t2+line[1], line[2]*t2+line[3], line[4]*t2+line[5]])
            return point1, point2
        elif discriminant < 0:
            raise ValueError('Discriminant less than 0.')

# Find the difference in time-of-arrival for the neutrinos and 2 detectors located at pos1 and pos2 and a supernova
# located at posSN
#
# Input: list/numpy array [x1, y1, z1], list/numpy array [x2, y2, z2], list/numpy array [x, y, z], number
#
# Output: number
#
def find_deltaT(pos1, pos2, posSN, c_light):
    local_posSN = copy.copy(posSN) # Define a local version of the SN position
    dT_max = distance_between_3d_points(pos1, pos2) / c_light
    relative_pos = pos2 - pos1
    relative_pos /= np.linalg.norm(relative_pos)
    local_posSN /= np.linalg.norm(local_posSN)
    return (dT_max * np.dot(local_posSN, relative_pos)) / (np.linalg.norm(local_posSN) * np.linalg.norm(relative_pos))

# Function that takes in quantities defining a circle in 3D space and returns a point on the circle at
# a user-given angle
#
# Input: list/numpy array [x1, y1, z1], list/numpy array [x2, y2, z2], number, number (angle in radians)
#
# Output: list/numpy array [x, y, z]
#
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
    return centre + radius * math.cos(theta) * a + radius * math.sin(theta) * b

# Function to produce a Skycoord object of long/lats of a 3D circle based on
# the ToA difference between 2 detectors
# The wrapped long/lats are returned
#
# Input: number, list/numpy array [x1, y1, z1], list/numpy array [x2, y2, z2], number, number, number
#
# Output: Skycoord array, Skycoord array
#
def circle_plot_points_astropy(N_points, pos1, pos2, c_light, dT, SN_dist):

    angle_list = np.linspace(0, 2 * np.pi, N_points)

    # Circle parameters
    theta = math.acos(dT * c_light / distance_between_3d_points(pos1, pos2))
    radius = SN_dist * math.sin(theta)
    circle_normal = pos2 - pos1
    circle_normal /= np.linalg.norm(circle_normal) # Normalize the plane normal vector
    circle_centre = SN_dist*math.cos(theta) * circle_normal
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
