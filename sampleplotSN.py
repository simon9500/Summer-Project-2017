import numpy as np
import math
import copy
import SNfunctions as SN
import sys
import random
from time import time
from ROOT import TCanvas, TH2D, gStyle
from scipy.optimize import curve_fit
from SNfunctions import distance_between_3d_points, point_on_3d_circle, lat_long_from_xyz

# This code samples from the SN neutrino distributions and produces a start time each .  The dt values are
# changed to account for this and the resulting circle points are added to a 2D aitoff histogram, which is saved
# at the end.

# For reference, the axes are defined as:
# x axis is along the line of 0 longitude and 0 latitude
# y axis is along the line of 90 longitude and 0 latitude
# z axis is along the line of 90 latitude (north pole)
# The origin is at the centre of earth

# Find a start time based on a time distribution (1ms bins) and the total number of events expected
# Sample a time from the distribution for each event and find the minimum time to an accuracy of 0.01ms
def find_start_time(total_events_scaled, total_events, time_val_list, diffrate_list, running_integral):
    start_time = 1000.  # Set start time to large number for the purpose of calculating a minimum
    for i in range(total_events_scaled):
        random_sample = random.random() * total_events
        first_encounter = 5000.
        for j in range(len(time_val_list)):
            if running_integral[j] >= random_sample:
                # Interpolate between the interval endpoints
                # Current integral is the integral up to the start of the interval
                # tvals and dr_interp are the interpolated data points over the 1ms interval
                if j == 0:
                    tvals = np.linspace(0., time_val_list[j], 101)
                    dr_interp = np.interp(tvals, [0., time_val_list[j]],
                                          [diffrate_list[j], diffrate_list[j]])
                    current_integral = 0.
                else:
                    tvals = np.linspace(time_val_list[j - 1], time_val_list[j], 101)
                    dr_interp = np.interp(tvals, [time_val_list[j - 1], time_val_list[j]],
                                          [diffrate_list[j - 1], diffrate_list[j]])
                    current_integral = running_integral[j - 1]
                # Run through integral and find the time value of the random sample to a precision of 0.01ms
                for n in range(1, len(tvals)):
                    current_integral += dr_interp[n] / 100000.
                    if current_integral >= random_sample:
                        first_encounter = tvals[n]
                        break
                break
        if first_encounter < start_time:
            start_time = first_encounter
    return start_time

# ALTERNATIVE METHOD - TAKES TOO LONG
# Find a start time based on a time distribution (0.01ms bins) and the total number of events expected
# Sample a time from the distribution for each event and find the minimum time to an accuracy of 0.01ms
def find_start_time_alternative_method(total_events_scaled, total_events,
                                       time_val_list, running_integral):
    start_time = 10000.  # Set start time to large number for the purpose of calculating a minimum
    for i in range(total_events_scaled):
        random_sample = random.random() * total_events
        first_encounter = 50.
        for j in range(len(time_val_list)):
            if running_integral[j] >= random_sample:
                first_encounter = time_val_list[j] # Round up the time value
                break
        if first_encounter < start_time:
            start_time = first_encounter
    return start_time

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

# Run the following code only if the code is being ran directly, ie. not being imported
if __name__ == '__main__':
    # Read names of output file and input file from command line
    if len(sys.argv) != 3:
        print "Wrong number of arguments."
        print "Usage: " + sys.argv[0] + " <(diffrate) input file> <(energy threshold) input file>"
        quit()
    else:
        infileName1 = sys.argv[1]
        infileName2 = sys.argv[2]

    beginning_time = time()

    # Open input file for reading
    infile1 = open(infileName1, "r")
    infile2 = open(infileName2, "r")

    # Initialise arrays
    time_val_list = []
    diffrate_list = []
    intrate_list = []
    enthreshold_list = []

    # Scaling factors for total events
    energy_threshold = 0.  # Energy threshold of the detector in keV
    detector_mass = 7.  # Detector mass in tonnes (7 for LZ)
    SN_dist_unit = 3.0857e+20  # 10kpc in m - all readings are relative to this distance and scale as 1/r^2
    SN_dist = 3.0857e+20  # 1kpc in m is 3.0857e+19
    SN_dist_event_scaling = 1. / (SN_dist / SN_dist_unit) ** 2  # Scaling factor for events by distance SN is from earth

    # Scaling factors for number of events from different detectors
    # We are assuming the differential rate has the same distribution for each detector, only the scaling is affected.
    # Hence, we can sample a time from the same distribution and only change the total events we expect using a
    # scaling factor.

    # For a SN at 8.5kpc:
    # LZ expects 261 events
    # SK expects 7000 events
    # SNO expects 300 events
    # GS (Icarus) expects 200 events
    # ICE CUBE is a different kind of detector so the events cannot be quantified in the same way
    scale_SK = 7000./261.
    scale_SNO = 300./261.
    scale_GS = 200./261.

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

    radius_earth = 6371000.  # Radius of earth in m
    c_light = 3.e+8  # Speed of light in m/s
    N_circle_points = 300  # Number of points to represent each circle

    # (x,y,z) coordinates of the detectors - found using website in Alex's sample.f code
    SK_pos = np.array([-3778227.70054802, 3500760.07816377, 3749077.00422199])
    SNO_pos = np.array([674206.505268229, -4348482.39235766, 4601529.68579595])
    LZ_pos = np.array([-1085869.32956348, -4437183.82092001, 4436246.60262743])
    GS_pos = np.array([4580523.90968328, 1105217.78652699, 4284183.43324446])
    ICE_pos = np.array([499.192610467315, -999.174518693858, -6356722.68768389])

    # Define SN latitude and longitude
    SN_lat = 45.
    SN_long = 45.

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
    dT_LZ_SK = dummy_dT_1  # ToA at LZ - ToA at SK
    dT_LZ_SNO = dummy_dT_2  # ToA at LZ - ToA at SNO
    dT_LZ_GS = dummy_dT_3  # ToA at LZ - ToA at GS
    dT_LZ_ICE = dummy_dT_4  # ToA at LZ - ToA at ICE
    dT_SK_ICE = dummy_dT_5  # ToA at SK - ToA at ICE
    dT_SK_SNO = dummy_dT_6  # ToA at SK - ToA at SNO
    dT_SK_GS = dummy_dT_7  # ToA at SK - ToA at GS
    dT_SNO_GS = dummy_dT_8  # ToA at SNO - ToA at GS
    dT_SNO_ICE = dummy_dT_9  # ToA at SNO - ToA at ICE
    dT_GS_ICE = dummy_dT_10  # ToA at GS - ToA at ICE

    # Find max delta t values between all detector pairs
    max_dT_LZ_SK = SN.distance_between_3d_points(LZ_pos, SK_pos) / c_light
    max_dT_LZ_SNO = SN.distance_between_3d_points(LZ_pos, SNO_pos) / c_light
    max_dT_LZ_GS = SN.distance_between_3d_points(LZ_pos, GS_pos) / c_light
    max_dT_LZ_ICE = SN.distance_between_3d_points(LZ_pos, ICE_pos) / c_light
    max_dT_SK_ICE = SN.distance_between_3d_points(SK_pos, ICE_pos) / c_light
    max_dT_SK_SNO = SN.distance_between_3d_points(SK_pos, SNO_pos) / c_light
    max_dT_SK_GS = SN.distance_between_3d_points(SK_pos, GS_pos) / c_light
    max_dT_SNO_GS = SN.distance_between_3d_points(SNO_pos, GS_pos) / c_light
    max_dT_SNO_ICE = SN.distance_between_3d_points(SNO_pos, ICE_pos) / c_light
    max_dT_GS_ICE = SN.distance_between_3d_points(GS_pos, ICE_pos) / c_light

    print "Input dt values (in ms)"
    print "LZ-SK = %f" % (1000 * dT_LZ_SK)
    print "LZ-SNO = %f" % (1000 * dT_LZ_SNO)
    print "LZ-GS = %f" % (1000 * dT_LZ_GS)
    print "LZ-ICE = %f" % (1000 * dT_LZ_ICE)
    print "SK-ICE = %f" % (1000 * dT_SK_ICE)
    print "SK-SNO = %f" % (1000 * dT_SK_SNO)
    print "SK-GS = %f" % (1000 * dT_SK_GS)
    print "SNO-GS = %f" % (1000 * dT_SNO_GS)
    print "SNO-ICE = %f" % (1000 * dT_SNO_ICE)
    print "GS-ICE = %f" % (1000 * dT_GS_ICE)

    canv = TCanvas("c", "c", 1200, 960)
    roothist = TH2D("ha", "Aitoff", 360, -180, 180, 180, -90, 90)

    # Fill arrays
    for line in infile1:
        tokens = line.split()
        time_val_list.append(float(tokens[0]))
        diffrate_list.append(float(tokens[1]))

    for line in infile2:
        tokens = line.split(",")
        enthreshold_list.append(float(tokens[0]))
        intrate_list.append(float(tokens[1]))

    # Close files
    infile1.close()
    infile2.close()

    # Convert lists to numpy arrays
    time_val_list = np.array(time_val_list)
    diffrate_list = np.array(diffrate_list)
    intrate_list = np.array(intrate_list)
    enthreshold_list = np.array(enthreshold_list)

    # Fit threshold data with an exponential
    def exp_func(x, a, b, c):
        return a * np.exp(-b * x) + c

    popt, pcov = curve_fit(exp_func, enthreshold_list, intrate_list)

    total_events = np.sum(diffrate_list) / 1000.  # Total integral of the file
    # Scale by SN distance, detector mass and the threshold energy of the detector
    # Scaling by mass is done by multiplying by the mass in tonnes of the detector
    # The total events scale off as 1/r^2, where r is the SN distance
    # The scaling by energy threshold is done by finding the value of the fitted function at the energy
    # threshold of the detector and using this as the baseline number of events expected for a 1T detector and 10kpc
    # SN.
    total_events_scaled = SN_dist_event_scaling * detector_mass * exp_func(energy_threshold, *popt)
    num_start_time = 100  # Number of start times to loop through
    num_events_LZ = 100  #int(total_events_scaled)
    num_events_SK = int(scale_SK * num_events_LZ)
    num_events_SNO = int(scale_SNO * num_events_LZ)
    num_events_GS = int(scale_GS * num_events_LZ)
    print "Number of start times sampled = %d" % num_start_time
    print "Total number of events = %d" % num_events_LZ
    print "%d events expected at LZ for a detector mass of %f and a SN distance of %fkpc" %(int(total_events_scaled),
                                                                                            detector_mass,
                                                                                            (SN_dist/SN_dist_unit)*10)
    outputfile_name = "sampleaitoffhist_" + str(num_events_LZ) + "LZevents_" + str(num_start_time) + "st.png"

    # Running integral (over the full 10s range) array - to save computation in the loop
    running_integral_1ms_bins = np.empty(len(diffrate_list))
    integral_sofar = 0.
    for m in range(len(diffrate_list)):
        integral_sofar += diffrate_list[m] / 1000.
        running_integral_1ms_bins[m] = integral_sofar

    """
    # Interpolate between each 1ms 100 times to produce values to bins of 0.01ms
    time_val_list_10micros_bins = np.linspace(1.00000005e-03, 10., len(time_val_list)*100)
    diffrate_list_10micros_bins = np.interp(time_val_list_10micros_bins, time_val_list, diffrate_list)
    running_integral_10micros_bins = np.empty(len(diffrate_list)*100)
    integral_sofar = 0.
    for m in range(len(running_integral_10micros_bins)):
        integral_sofar += diffrate_list_10micros_bins[m] / 100000.
        running_integral_10micros_bins = integral_sofar
    """

    # Main loop
    # Produce circles according to the adjusted start times and add them to an Aitoff histogram
    for k in range(num_start_time):
        start_time_LZ = find_start_time(num_events_LZ, total_events, time_val_list,
                                        diffrate_list, running_integral_1ms_bins)
        start_time_SK = find_start_time(num_events_SK, total_events, time_val_list,
                                        diffrate_list, running_integral_1ms_bins)
        start_time_SNO = find_start_time(num_events_SNO, total_events, time_val_list,
                                         diffrate_list, running_integral_1ms_bins)
        start_time_GS = find_start_time(num_events_GS, total_events, time_val_list,
                                        diffrate_list, running_integral_1ms_bins)
        """
        print "LZ Start time = %fms." %(start_time_LZ * 1000)
        print "SNO Start time = %fms." % (start_time_SNO * 1000)
        print "SK Start time = %fms." % (start_time_SK * 1000)
        print "GS Start time = %fms." % (start_time_GS * 1000)
        """

        # Apply start time uncertainty in LZ to dT values
        dT1 = dT_LZ_GS + start_time_LZ - start_time_GS
        dT2 = dT_LZ_SNO + start_time_LZ - start_time_SNO
        dT3 = dT_LZ_SK + start_time_LZ - start_time_SK
        dT4 = dT_LZ_ICE + start_time_LZ
        dT5 = dT_SK_SNO + start_time_SK - start_time_SNO
        dT6 = dT_SK_GS + start_time_SK - start_time_GS
        dT7 = dT_SK_ICE + start_time_SK
        dT8 = dT_SNO_GS + start_time_SNO - start_time_GS
        dT9 = dT_SNO_ICE + start_time_SNO
        dT10 = dT_GS_ICE + start_time_GS

        # Create set of circle points (lat/long) using adjusted dT values
        # The conditional statements are included to ensure the adjusted dT vales are appropriate
        if abs(dT1) <= max_dT_LZ_GS:
            lat_LZ_GS, long_LZ_GS = circle_plot_points(N_circle_points,
                                                       LZ_pos, GS_pos, c_light, dT1, SN_dist)
            # Fill histogram with circle long/lats (weight = 1)
            for i in range(N_circle_points):
                roothist.Fill(long_LZ_GS[i], lat_LZ_GS[i], 1)
        if abs(dT2) <= max_dT_LZ_SNO:
            lat_LZ_SNO, long_LZ_SNO = circle_plot_points(N_circle_points,
                                                         LZ_pos, SNO_pos, c_light, dT2, SN_dist)
            # Fill histogram with circle long/lats (weight = 1)
            for i in range(N_circle_points):
                roothist.Fill(long_LZ_SNO[i], lat_LZ_SNO[i], 1)
        if abs(dT3) <= max_dT_LZ_SK:
            lat_LZ_SK, long_LZ_SK = circle_plot_points(N_circle_points,
                                                       LZ_pos, SK_pos, c_light, dT3, SN_dist)
            # Fill histogram with circle long/lats (weight = 1)
            for i in range(N_circle_points):
                roothist.Fill(long_LZ_SK[i], lat_LZ_SK[i], 1)
        if abs(dT4) <= max_dT_LZ_ICE:
            lat_LZ_ICE, long_LZ_ICE = circle_plot_points(N_circle_points,
                                                         LZ_pos, ICE_pos, c_light, dT4, SN_dist)
            # Fill histogram with circle long/lats (weight = 1)
            for i in range(N_circle_points):
                roothist.Fill(long_LZ_ICE[i], lat_LZ_ICE[i], 1)
        if abs(dT5) <= max_dT_SK_SNO:
            lat_SK_SNO, long_SK_SNO = circle_plot_points(N_circle_points,
                                                         SK_pos, SNO_pos, c_light, dT5, SN_dist)
            # Fill histogram with circle long/lats (weight = 1)
            for i in range(N_circle_points):
                roothist.Fill(long_SK_SNO[i], lat_SK_SNO[i], 1)
        if abs(dT6) <= max_dT_SK_GS:
            lat_SK_GS, long_SK_GS = circle_plot_points(N_circle_points,
                                                       SK_pos, GS_pos, c_light, dT6, SN_dist)
            # Fill histogram with circle long/lats (weight = 1)
            for i in range(N_circle_points):
                roothist.Fill(long_SK_GS[i], lat_SK_GS[i], 1)
        if abs(dT7) <= max_dT_SK_ICE:
            lat_SK_ICE, long_SK_ICE = circle_plot_points(N_circle_points,
                                                         SK_pos, ICE_pos, c_light, dT7, SN_dist)
            # Fill histogram with circle long/lats (weight = 1)
            for i in range(N_circle_points):
                roothist.Fill(long_SK_ICE[i], lat_SK_ICE[i], 1)
        if abs(dT8) <= max_dT_SNO_GS:
            lat_SNO_GS, long_SNO_GS = circle_plot_points(N_circle_points,
                                                         SNO_pos, GS_pos, c_light, dT8, SN_dist)
            # Fill histogram with circle long/lats (weight = 1)
            for i in range(N_circle_points):
                roothist.Fill(long_SNO_GS[i], lat_SNO_GS[i], 1)
        if abs(dT9) <= max_dT_SNO_ICE:
            lat_SNO_ICE, long_SNO_ICE = circle_plot_points(N_circle_points,
                                                           SNO_pos, ICE_pos, c_light, dT9, SN_dist)
            # Fill histogram with circle long/lats (weight = 1)
            for i in range(N_circle_points):
                roothist.Fill(long_SNO_ICE[i], lat_SNO_ICE[i], 1)
        if abs(dT10) <= max_dT_GS_ICE:
            lat_GS_ICE, long_GS_ICE = circle_plot_points(N_circle_points,
                                                         GS_pos, ICE_pos, c_light, dT10, SN_dist)
            # Fill histogram with circle long/lats (weight = 1)
            for i in range(N_circle_points):
                roothist.Fill(long_GS_ICE[i], lat_GS_ICE[i], 1)

        print "%d" % (k + 1)

    gStyle.SetOptStat(0)  # Remove stat boxes
    roothist.Draw("aitoff")
    canv.SaveAs(outputfile_name)

    end_time = time()

    print "Time to run program = %fs" %(end_time - beginning_time)

    raise SystemExit
