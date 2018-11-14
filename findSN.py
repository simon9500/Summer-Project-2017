import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import SNfunctions as SN

# This code finds the location in the sky of a SN based on the neutrino ToA differences observed in different
# detectors, then plots the circles on an Aitoff projection.

# For reference, the axes are defined as:
# x axis is along the line of 0 longitude and 0 latitude
# y axis is along the line of 90 longitude and 0 latitude
# z axis is along the line of 90 latitude (north pole)
# The origin is at the centre of earth

# Run the following code only if the code is being ran directly, ie. not being imported
if __name__=='__main__':
    # Latitudes are positive for N, and negative for S.  Longitudes are positive for E and negative for W.
    # Changed to agree with Alex's values
    SK_lat = 36.2333324
    SK_long = 137.1833326
    SNO_lat = 46.4719
    SNO_long = -81.1868
    LZ_lat = 44.3533
    LZ_long = -103.7511861
    GS_lat = 42.4691
    GS_long = 13.5654
    ICE_lat = -89.99
    ICE_long = -63.453056

    SN_dist = 1.e+14  # Distance SN is from earth - 10kpc is 3.0857e+20 in m
    radius_earth = 6371000.  # Radius of earth in m
    c_light = 3.e+8  # Speed of light in m/s
    N_circle_points = 500 # Number of points to represent each circle

    # (x,y,z) coordinates of the detectors - found using website in Alex's sample.f code
    SK_pos = np.array([-3778227.70054802,3500760.07816377,3749077.00422199])
    SNO_pos = np.array([674206.505268229,-4348482.39235766,4601529.68579595])
    LZ_pos = np.array([-1085869.32956348,-4437183.82092001,4436246.60262743])
    GS_pos = np.array([4580523.90968328,1105217.78652699,4284183.43324446])
    ICE_pos = np.array([499.192610467315,-999.174518693858,-6356722.68768389])

    # Define SN latitude and longitude
    SN_lat = -20.
    SN_long = 155.

    # Choose a SN position
    posSN = SN.find_xyz_from_latlong(SN_lat, SN_long, SN_dist)
    dummy_dT_1 = SN.find_deltaT(LZ_pos, SK_pos, posSN, c_light)
    dummy_dT_2 = SN.find_deltaT(LZ_pos, SNO_pos, posSN, c_light)
    dummy_dT_3 = SN.find_deltaT(LZ_pos, GS_pos, posSN, c_light)
    dummy_dT_4 = SN.find_deltaT(LZ_pos, ICE_pos, posSN, c_light)
    dummy_dT_5 = SN.find_deltaT(SK_pos, ICE_pos, posSN, c_light)
    dummy_dT_6 = SN.find_deltaT(SK_pos, SNO_pos, posSN, c_light)
    dummy_dT_7 = SN.find_deltaT(SK_pos, GS_pos, posSN, c_light)
    dummy_dT_8 = SN.find_deltaT(SNO_pos, GS_pos, posSN, c_light)
    dummy_dT_9 = SN.find_deltaT(SNO_pos, ICE_pos, posSN, c_light)
    dummy_dT_10 = SN.find_deltaT(GS_pos, ICE_pos, posSN, c_light)

    # Define the delta t values for the detectors
    dT_LZ_SK = dummy_dT_1 # ToA at LZ - ToA at SK
    dT_LZ_SNO = dummy_dT_2 # ToA at LZ - ToA at SNO
    dT_LZ_GS = dummy_dT_3 # ToA at LZ - ToA at GS
    dT_LZ_ICE = dummy_dT_4 # ToA at LZ - ToA at ICE
    dT_SK_ICE = dummy_dT_5 # ToA at SK - ToA at ICE
    dT_SK_SNO = dummy_dT_6 # ToA at SK - ToA at SNO
    dT_SK_GS = dummy_dT_7 # ToA at SK - ToA at GS
    dT_SNO_GS = dummy_dT_8 # ToA at SNO - ToA at GS
    dT_SNO_ICE = dummy_dT_9 # ToA at SNO - ToA at ICE
    dT_GS_ICE = dummy_dT_10 # ToA at GS - ToA at ICE
    print "------------------------------------------------------"
    print "Input dt values (in ms)"
    print "LZ-SK = %f" %(1000 * dT_LZ_SK)
    print "LZ-SNO = %f" % (1000 * dT_LZ_SNO)
    print "LZ-GS = %f" % (1000 * dT_LZ_GS)
    print "LZ-ICE = %f" % (1000 * dT_LZ_ICE)
    print "SK-ICE = %f" % (1000 * dT_SK_ICE)
    print "SK-SNO = %f" % (1000 * dT_SK_SNO)
    print "SK-GS = %f" % (1000 * dT_SK_GS)
    print "SNO-GS = %f" % (1000 * dT_SNO_GS)
    print "SNO-ICE = %f" % (1000 * dT_SNO_ICE)
    print "GS-ICE = %f" % (1000 * dT_GS_ICE)

    # Produce plane equations for all circles created from the detector pairs
    plane_1 = SN.produce_plane(LZ_pos, SK_pos, dT_LZ_SK, SN_dist, c_light)
    plane_2 = SN.produce_plane(LZ_pos, SNO_pos, dT_LZ_SNO, SN_dist, c_light)
    plane_3 = SN.produce_plane(LZ_pos, GS_pos, dT_LZ_GS, SN_dist, c_light)
    plane_4 = SN.produce_plane(LZ_pos, ICE_pos, dT_LZ_ICE, SN_dist, c_light)
    plane_5 = SN.produce_plane(SK_pos, ICE_pos, dT_SK_ICE, SN_dist, c_light)
    plane_6 = SN.produce_plane(SK_pos, SNO_pos, dT_SK_SNO, SN_dist, c_light)
    plane_7 = SN.produce_plane(SK_pos, GS_pos, dT_SK_GS, SN_dist, c_light)
    plane_8 = SN.produce_plane(SNO_pos, GS_pos, dT_SNO_GS, SN_dist, c_light)
    plane_9 = SN.produce_plane(SNO_pos, ICE_pos, dT_SNO_ICE, SN_dist, c_light)
    plane_10 = SN.produce_plane(GS_pos, ICE_pos, dT_GS_ICE, SN_dist, c_light)

    # Define equation of sphere centred on earth
    sphere = np.array([0., 0., 0., SN_dist])

    """
    Finally, we calculate the position of the SN.
    """

    plane_list = [plane_1, plane_2, plane_3, plane_4, plane_5, plane_6, plane_7, plane_8, plane_9, plane_10]
    intersection_points_list = []

    for i in range(len(plane_list)):
        for j in range(len(plane_list)):
            if j > i:
                current_line_int = SN.plane_intersection(plane_list[i], plane_list[j])
                point1, point2 = SN.intersection_2planes_sphere(current_line_int, sphere)
                points = np.array([SN.lat_long_from_xyz(point1), SN.lat_long_from_xyz(point2)])
                intersection_points_list.append(points)

    intersection_points_list = np.array(intersection_points_list)
    print len(intersection_points_list)
    print "------------------------------------------------------"
    for pts in intersection_points_list:
        for p in pts:
            print "Latitude = %f, Longitude = %f" %(p[0], p[1])
    print "------------------------------------------------------"


    # Finding the line of intersection between a selection of planes
    intersect_line_1 = SN.plane_intersection(plane_1, plane_2)
    intersect_line_2 = SN.plane_intersection(plane_2, plane_3)
    intersect_line_3 = SN.plane_intersection(plane_1, plane_3)
    intersect_line_4 = SN.plane_intersection(plane_1, plane_4)
    intersect_line_5 = SN.plane_intersection(plane_1, plane_5)
    intersect_line_6 = SN.plane_intersection(plane_2, plane_6)
    intersect_line_7 = SN.plane_intersection(plane_3, plane_7)
    intersect_line_8 = SN.plane_intersection(plane_9, plane_8)
    intersect_line_9 = SN.plane_intersection(plane_10, plane_4)
    intersect_line_10 = SN.plane_intersection(plane_7, plane_5)

    # Finding each pair of circle-circle intersection points
    point1, point2 = SN.intersection_2planes_sphere(intersect_line_1, sphere)
    point3, point4 = SN.intersection_2planes_sphere(intersect_line_2, sphere)
    point5, point6 = SN.intersection_2planes_sphere(intersect_line_3, sphere)
    point7, point8 = SN.intersection_2planes_sphere(intersect_line_4, sphere)
    point9, point10 = SN.intersection_2planes_sphere(intersect_line_5, sphere)
    point11, point12 = SN.intersection_2planes_sphere(intersect_line_6, sphere)
    point13, point14 = SN.intersection_2planes_sphere(intersect_line_7, sphere)
    point15, point16 = SN.intersection_2planes_sphere(intersect_line_8, sphere)
    point17, point18 = SN.intersection_2planes_sphere(intersect_line_9, sphere)
    point19, point20 = SN.intersection_2planes_sphere(intersect_line_10, sphere)

    # Arrays to hold the circle-circle pair intersection locations (in latitude and longitude)
    points1 = np.array([SN.lat_long_from_xyz(point1), SN.lat_long_from_xyz(point2)])
    points2 = np.array([SN.lat_long_from_xyz(point3), SN.lat_long_from_xyz(point4)])
    points3 = np.array([SN.lat_long_from_xyz(point5), SN.lat_long_from_xyz(point6)])
    points4 = np.array([SN.lat_long_from_xyz(point7), SN.lat_long_from_xyz(point8)])
    points5 = np.array([SN.lat_long_from_xyz(point9), SN.lat_long_from_xyz(point10)])
    points6 = np.array([SN.lat_long_from_xyz(point11), SN.lat_long_from_xyz(point12)])
    points7 = np.array([SN.lat_long_from_xyz(point13), SN.lat_long_from_xyz(point14)])
    points8 = np.array([SN.lat_long_from_xyz(point15), SN.lat_long_from_xyz(point16)])
    points9 = np.array([SN.lat_long_from_xyz(point17), SN.lat_long_from_xyz(point18)])
    points10 = np.array([SN.lat_long_from_xyz(point19), SN.lat_long_from_xyz(point20)])

    # Initialise variable to hold SN latitude/longitude
    SN_position = [0., 0.]

    # Locate the SN
    for p1 in points1:
        for p2 in points2:
            for p3 in points3:
                for p4 in points4:
                    for p5 in points5:
                        for p6 in points6:
                            for p7 in points7:
                                for p8 in points8:
                                    for p9 in points9:
                                        for p10 in points10:
                                            sep_lat_1 = abs(p1[0] - p2[0])
                                            sep_lat_2 = abs(p1[0] - p3[0])
                                            sep_lat_3 = abs(p1[0] - p4[0])
                                            sep_lat_4 = abs(p1[0] - p5[0])
                                            sep_lat_5 = abs(p1[0] - p6[0])
                                            sep_lat_6 = abs(p1[0] - p7[0])
                                            sep_lat_7 = abs(p1[0] - p8[0])
                                            sep_lat_8 = abs(p1[0] - p9[0])
                                            sep_lat_9 = abs(p1[0] - p10[0])
                                            sep_long_1 = abs(p1[1] - p2[1])
                                            sep_long_2 = abs(p1[1] - p3[1])
                                            sep_long_3 = abs(p1[1] - p4[1])
                                            sep_long_4 = abs(p1[1] - p5[1])
                                            sep_long_5 = abs(p1[1] - p6[1])
                                            sep_long_6 = abs(p1[1] - p7[1])
                                            sep_long_7 = abs(p1[1] - p8[1])
                                            sep_long_8 = abs(p1[1] - p9[1])
                                            sep_long_9 = abs(p1[1] - p10[1])
                                            if sep_lat_1<0.5 and sep_lat_2<0.5 and \
                                                    sep_long_1<0.5 and sep_long_2<0.5 and \
                                                    sep_lat_3<0.5 and sep_long_3<0.5 and \
                                                    sep_lat_4<0.5 and sep_long_4<0.5 and \
                                                    sep_lat_5<0.5 and sep_long_5<0.5 and \
                                                    sep_lat_6 < 0.5 and sep_long_6 < 0.5 and \
                                                    sep_lat_7 < 0.5 and sep_long_7 < 0.5 and \
                                                    sep_lat_8 < 0.5 and sep_long_8 < 0.5 and \
                                                    sep_lat_9 < 0.5 and sep_long_9 < 0.5:
                                                SN_position = [(p1[0]+p2[0]+p3[0]+p4[0]+p5[0]+
                                                                p6[0]+p7[0]+p8[0]+p9[0]+p10[0])/10.,
                                                               (p1[1]+p2[1]+p3[1]+p4[1]+p5[1]+
                                                                p6[1]+p7[1]+p8[1]+p9[1]+p10[1])/10.]
                                            elif p1[0] == p2[0] == p3[0] == p4[0] == p5[0] == p6[0] == p7[0] == p8[0] \
                                                 == p9[0] == p10[0] == 90.:
                                                SN_position = [90.,0.]
                                            elif p1[0] == p2[0] == p3[0] == p4[0] == p5[0] == p6[0] == p7[0] == p8[0] \
                                                 == p9[0] == p10[0] == -90.:
                                                SN_position = [-90.,0.]

    print "The SN is on latitude %f and longitude %f" %(SN_position[0], SN_position[1])
    print "------------------------------------------------------"

    """
    Plot the circles from each detector pair, the detectors and the SN.
    """

    c_Det = SkyCoord(ra=[LZ_long, GS_long, SK_long, SNO_long, ICE_long]*u.degree,
                     dec=[LZ_lat, GS_lat, SK_lat, SNO_lat, ICE_lat]*u.degree
                    , frame='icrs')
    c_SN = SkyCoord(ra=[SN_position[1]]*u.degree, dec=[SN_position[0]]*u.degree, frame='icrs')
    
    ra_rad_Det = c_Det.ra.wrap_at(180 * u.deg).radian
    dec_rad_Det = c_Det.dec.radian
    ra_rad_SN = c_SN.ra.wrap_at(180 * u.deg).radian
    dec_rad_SN = c_SN.dec.radian
    
    
    ra_rad_LZ_GS, dec_rad_LZ_GS = SN.circle_plot_points_astropy(N_circle_points,
                                                             LZ_pos, GS_pos, c_light, dT_LZ_GS, SN_dist)
    ra_rad_LZ_SNO, dec_rad_LZ_SNO = SN.circle_plot_points_astropy(N_circle_points,
                                                               LZ_pos, SNO_pos, c_light, dT_LZ_SNO, SN_dist)
    ra_rad_LZ_SK, dec_rad_LZ_SK = SN.circle_plot_points_astropy(N_circle_points,
                                                               LZ_pos, SK_pos, c_light, dT_LZ_SK, SN_dist)
    ra_rad_LZ_ICE, dec_rad_LZ_ICE = SN.circle_plot_points_astropy(N_circle_points,
                                                               LZ_pos, ICE_pos, c_light, dT_LZ_ICE, SN_dist)
    ra_rad_SK_SNO, dec_rad_SK_SNO = SN.circle_plot_points_astropy(N_circle_points,
                                                               SK_pos, SNO_pos, c_light, dT_SK_SNO, SN_dist)
    ra_rad_SK_GS, dec_rad_SK_GS = SN.circle_plot_points_astropy(N_circle_points,
                                                               SK_pos, GS_pos, c_light, dT_SK_GS, SN_dist)
    ra_rad_SK_ICE, dec_rad_SK_ICE = SN.circle_plot_points_astropy(N_circle_points,
                                                               SK_pos, ICE_pos, c_light, dT_SK_ICE, SN_dist)
    ra_rad_SNO_GS, dec_rad_SNO_GS = SN.circle_plot_points_astropy(N_circle_points,
                                                               SNO_pos, GS_pos, c_light, dT_SNO_GS, SN_dist)
    ra_rad_SNO_ICE, dec_rad_SNO_ICE = SN.circle_plot_points_astropy(N_circle_points,
                                                               SNO_pos, ICE_pos, c_light, dT_SNO_ICE, SN_dist)
    ra_rad_GS_ICE, dec_rad_GS_ICE = SN.circle_plot_points_astropy(N_circle_points,
                                                               GS_pos, ICE_pos, c_light, dT_GS_ICE, SN_dist)

    # As a last step we set up the plotting environment with matplotlib using the
    # Aitoff projection with a specific title, a grid, filled circles as markers with
    # a markersize of 2 and an alpha value of 0.3.
    plt.figure(figsize=(8,4.2))
    plt.subplot(111, projection="aitoff")
    plt.title("Aitoff projection", y=1.08)
    plt.grid(True)

    # Plot all circles
    plt.plot(ra_rad_LZ_GS, dec_rad_LZ_GS, 'o', color='green', markersize=2., alpha=0.6, label='LZ-GS circle')
    plt.plot(ra_rad_LZ_SK, dec_rad_LZ_SK, 'o', color='red', markersize=2., alpha=0.6, label ='LZ-SK circle')
    plt.plot(ra_rad_LZ_SNO, dec_rad_LZ_SNO, 'o', color='blue', markersize=2., alpha=0.6, label ='LZ-SNO circle')
    plt.plot(ra_rad_LZ_ICE, dec_rad_LZ_ICE, 'o', color='yellow', markersize=2., alpha=0.6, label='LZ-ICE circle')
    plt.plot(ra_rad_SK_SNO, dec_rad_SK_SNO, 'o', color='cyan', markersize=2., alpha=0.6, label='SK-SNO circle')
    plt.plot(ra_rad_SK_GS, dec_rad_SK_GS, 'o', color='magenta', markersize=2., alpha=0.6, label='SK-GS circle')
    plt.plot(ra_rad_SK_ICE, dec_rad_SK_ICE, 'o', color='#C5AD8E', markersize=2., alpha=0.6, label='SK-ICE circle')
    plt.plot(ra_rad_SNO_GS, dec_rad_SNO_GS, 'o', color='#35E781', markersize=2., alpha=0.6, label='SNO-GS circle')
    plt.plot(ra_rad_SNO_ICE, dec_rad_SNO_ICE, 'o', color='#C31AE5', markersize=2., alpha=0.6, label='SNO-ICE circle')
    plt.plot(ra_rad_GS_ICE, dec_rad_GS_ICE, 'o', color='#A0F40D', markersize=2., alpha=0.6, label='GS-ICE circle')

    # Plot all detectors
    plt.plot(ra_rad_Det, dec_rad_Det, 'o', color='black', markersize=4, alpha=1.0, label='Detectors')

    # Plot the SN
    plt.plot(ra_rad_SN, dec_rad_SN, 'o', color='orange', markersize=6, alpha=1.0, label='SN')
    plt.subplots_adjust(top=0.95, bottom=0.0)
    plt.legend(framealpha=0.3, bbox_to_anchor=(0.87, 0.75))
    plt.savefig('sky_example.png', bbox_inches='tight')
    plt.show()

    raise SystemExit
