'''
TODO
[x] Add 3D
[x] Add time_slice identification to each photon
[x] Examine data and identify source of abnormally large vectors.
[x] Examine why the direction seems to not be changing in a compton scattering
[x] Add polarisation
[x] Complete compton process
    [x] Energy change
    [x] Polarisation change
[x] Detection process
    [x] Where detection
        [x] Make faster
    [x] Blurring
        [x] Energy blurring
            [x] Fix error
        [x] Time blurring
            [x] Add as input
    [x] Store detection
    [x] Store all photons
        [x] check memory
    [x] Check why no detection is happening
[x] Pairing and noting event
    [x] Energy filter
    [x] Polarisation pairing
    [x] Test and debug pairing
[x] See if looping over time windows is necessary at all 
    I think it is, because when we do the actual simulation we'll have too many photons overall to complete in one go, memory constraints.
[x] Allow input for variables
[x] Check for true counts
[x] Vary source activity
[x] Plot fraction of correct coincidences vs activity

[] Add comments
[] Try implementing actual experiment method (pulse detection)
[] Root file for data
[x] Implement exercises 2 and 3
    [x] Different source location
        [x] Does compton need to be changed?
        [x] Make plots
    [x] Absorber
        [x] Absorb photons
        [x] Add exist checks to compton
        [x] " to detection
        [x] Add input options for parameters
        [Nope] Add varying absorption radius
'''


import numpy as np
# Imports a custom vector library
# The file vector_custom.py needs to be in the same directory or in the path.
import vector_custom as vec
import scipy as sc
import matplotlib.pyplot as plt
import copy
import time
import math
import random
from mpl_toolkits.mplot3d import axes3d, Axes3D

# Constants
hbar = 4.135667696e-15 # eV / Hz
# from scipy.constants import hbar
c = 299792458e3 # mm/s
# Not really necessary since we normalise the probability distribution so it cancels out
classical_electron_radius = 1 # 2.8179403262e-15 # m

# To take inputs or not
take_inputs = True
input_times_activities = False
input_compton_parameters = False
input_absorber_parameters = False
input_detector_parameters = False
# Make interactive plot?
make_plots = True

if __name__=="__main__" and take_inputs:
    if input("Would you like to use default values for time window, source activity? (y/n)")=='n':
        input_times_activities = True
    if input("Would you like to use default parameters for the compton ring? (y/n)")=='n':
        input_compton_parameters = True
    if input("Would you like to use default parameters for the absorbing material? (y/n)")=='n':
        input_absorber_parameters = True
    if input("Would you like to use default parameters for the detector? (y/n)")=='n':
        input_detector_parameters = True
    if input("Would you like to draw plots? (y/n)")=='y':
        make_plots = True

# The code is written so that it could be imported into another file and used as well.
# Hence there are functions available to set all input variables.
def define_make_plots(x:bool):
    # Need global keyword to make changes to global variables
    global make_plots
    make_plots = x

# Times and source activities
time_unit = 1 # nanoseconds
# Simulate one time window at a time because of memory constraints
time_window = 10 # nanoseconds
total_time = 10**4 # nanoseconds
source_activity_slice = 0.8 # nanoBecquerel - disintegrations per nanosecond

if input_times_activities:
    time_window = float(input("Coincidence time window (in nanoseconds) : "))
    total_time = float(input("Total time of simulation (in seconds) : "))*10**9 # To get in nanoseconds
n_windows = total_time//time_window

if input_times_activities:
    source_activity_slice = float(input("Source activity (in disintegrations per nanosecond) : "))
source_activity_window = source_activity_slice*time_window/time_unit # source activity per time window

def define_times_and_source_activities(time_unit_input, time_window_input, total_time_input, source_activity_slice_input):
    global time_unit,time_window, total_time, source_activity_slice, n_windows, source_activity_window
    time_unit = time_unit_input # nanoseconds

    time_window = time_window_input # nanoseconds
    total_time = total_time_input # nanoseconds
    source_activity_slice = source_activity_slice_input # nanoBecquerel - disintegrations per nanosecond

    n_windows = total_time//time_window
    source_activity_window = source_activity_slice*time_window/time_unit


# Compton ring
compton_probability = 0.2
compton_radius = 100 # mm
compton_width = 10 # mm

if input_compton_parameters:
    # Input probability and check for validity
    while True:
        compton_probability = float(input("Probability of a compton scattering : "))
        if 0<=compton_probability<=1:
            break
        print('Please enter a valid probability, between 0 and 1')
    compton_radius = float(input("Compton ring radius (in cm)"))*10 # mm
    compton_width = float(input("Compton ring width (in mm)")) # mm

def define_compton_parameters(compton_probability_input,compton_radius_input,compton_width_input):
    global compton_probability,compton_radius,compton_width
    compton_probability = compton_probability_input
    compton_radius = compton_radius_input
    compton_width = compton_width_input

# Absorbing material in compton ring
attenuation = 0.0096 # 0.0096 for water # absorbtions per mm
absorber_radius = compton_radius
absorber_width = compton_width

if input_absorber_parameters:
    attenuation = float(input("Attenuation in absorbtions per cm : "))/10 # absorbtions per mm
    if input("Use compton radius for absorber radius? (y/n)")=='n':
        absorber_radius = float(input("Absorbing disc radius (in cm)"))*10 # mm
    if input("Use compton width for absorber width? (y/n)")=='n':
        absorber_width = float(input("Absorbing disc width (in mm)")) # mm

def define_absorber_parameters(attenuation_input,
                               use_compton_for_absorber_radius:bool=True,
                               use_compton_for_absorber_width:bool=True,
                               absorber_radius_input=compton_radius,
                               absorber_width_input=compton_width):
    global attenuation,absorber_radius,absorber_width
    attenuation = attenuation_input
    if not use_compton_for_absorber_radius:
        absorber_radius = absorber_radius_input
    if not use_compton_for_absorber_width:
        absorber_width = absorber_width_input

# Detector
detector_radius = 200 # mm
detector_width = 10 # mm
detector_energy_cutoff_low = 200 # KeV
detector_energy_cutoff_high = 900 # KeV
timing_resolution = 1 # nanosec # FWHM of time blurring

if input_detector_parameters:
    detector_radius = float(input("Detector ring radius (in cm)"))*10 # mm
    detector_width = float(input("Detector's height along it's axis (in mm)")) # mm
    print("Detector's energy cutoffs (in KeV) : ")
    detector_energy_cutoff_low = float(input("Lower cutoff"))
    detector_energy_cutoff_high = float(input("Higher cutoff"))
    timing_resolution = float(input("Time Blurring FWHM / Timing resolution (in ns) : "))

def define_detector_parameters(detector_radius_input, detector_width_input, detector_energy_cutoff_low_input, detector_energy_cutoff_high_input, timing_resolution_input):
    global detector_radius,detector_width,detector_energy_cutoff_low,detector_energy_cutoff_high,timing_resolution
    detector_radius = detector_radius_input
    detector_width = detector_width_input
    detector_energy_cutoff_low = detector_energy_cutoff_low_input
    detector_energy_cutoff_high = detector_energy_cutoff_high_input
    timing_resolution = timing_resolution_input

# Unit vectors
# All vectors are lists in this program
hatz = [0,0,1]
hatx = [1,0,0]

def intersection(start_loc, direction, radius = detector_radius, shape = 0, **kwargs):
    '''
    Find where a ray with a given starting location and direction intersects a surface (detector/scatterer).
    Shape determines ring-shaped or spherical surface. Other shapes can be added later.
    Returns intersection point or None if no intersection.
    '''
    
    if shape==0: # A ring-shaped detector/scatterer with a fixed radius
        # Width only needed ring, so entered as a keyword argument (kwarg)
        try:
            width = kwargs['width']
        except KeyError:
            width = detector_width 

        # Want to solve the quadratic (start_loc[0]+direction[0]*t)**2 + (start_loc[1]+direction[1]*t)**2 - detector_radius**2 = 0
        # t is a parameter measuring how far along direction do we move from start_loc
        coef_a = direction[0]**2 + direction[1]**2
        coef_b = 2*(start_loc[0]*direction[0]+start_loc[1]*direction[1])
        coef_c = start_loc[0]**2 + start_loc[1]**2 - radius**2
        sqrt_determinant = np.sqrt(coef_b**2 - 4*coef_a*coef_c)
        # Quadratic formula
        sols = (-coef_b + np.array([sqrt_determinant,-sqrt_determinant]))/2/coef_a
        # Filter for moving along direction and not opposite to it
        solution = list(filter((lambda x: x>=0),sols))[0]


        # Also need to check if the z-coordinate is small enough that an intersection happens
        if abs(start_loc[2]+direction[2]*solution) <= width/2:
            return vec.add(start_loc,np.array(direction)*solution)
        return None
    
    if shape==1: # A spherical detector/scatterer with a fixed radius
        # Want to solve the quadratic sum_i ((start_loc[i]+direction[i]*t)**2) - detector_radius**2 = 0
        coef_a = direction[0]**2 + direction[1]**2 + direction[2]**2
        coef_b = 2*(start_loc[0]*direction[0]+start_loc[1]*direction[1] + start_loc[2]*direction[2])
        coef_c = start_loc[0]**2 + start_loc[1]**2 + start_loc[2]**2 - radius**2
        sqrt_determinant = np.sqrt(coef_b**2 - 4*coef_a*coef_c)

        sols = (-coef_b + np.array([sqrt_determinant,-sqrt_determinant]))/2/coef_a
        solution = list(filter((lambda x: x>=0),sols))[0]

        return vec.add(start_loc,np.array(direction)*solution)


# Without polarisation case - incomplete
# def Klein_Nishina_Distribution(theta, relative_energy):
#     if abs(theta)>np.pi:
#         return 0
#     elif abs(theta)<=1/np.sqrt(relative_energy):
#         return classical_electron_radius**2
#     else:
#         return

# With Polarisation, Klein-Nishima distribution
def Differential_Cross_Section(theta, phi, relative_energy):
    '''
    Calculate the differential cross section for compton scattering at a given energy and scattering angles.
    Theta is the scattering angle with respect to direction of motion.
    Phi is the scattering angle with respect to the polarisation vector in the plane normal to original direction of motion.
    Relative energy is energy per 511 KeV.
    '''
    # Phi has to be relative to polarisation angle
    lambda_ratio = 1/(1+relative_energy*(1-np.cos(theta)))
    return (classical_electron_radius*lambda_ratio)**2/2*(lambda_ratio+1/lambda_ratio-2*np.sin(theta)**2*np.cos(phi)**2)
    
def Sampling_Klein_Nishima(relative_energy):
    '''
    This function returns scattering angles theta, phi with probabilities given by the Klein-Nishima distribution.
    It uses rejection sampling, a method used to sample arbitrary probability distributions.
    Relative energy is energy per 511 KeV.
    '''
    # Rejection sampling:
    # First, choose theta based on a Gaussian and phi uniformly
    # Do probability check with differential cross section (times dOmega) and repeat if not true

    # Total cross section for normalisation
    if relative_energy==1.0: # The photons which aren't compton scattered will have exactly 1 relative energy
        total_cross_section = 3.608457130285029 # Precaculated for the common case
    else:
        total_cross_section = sc.integrate.dblquad(lambda theta, phi: Differential_Cross_Section(theta, phi, relative_energy)*np.sin(theta),
                                                   0,2*np.pi,
                                                   lambda x:0,lambda x:np.pi)[0]
    
    # The Klein-Nishima probability distribution obtained by normalising differential cross section * infinitesimal solid angle
    KN_PDF = lambda theta, phi : Differential_Cross_Section(theta, phi, relative_energy)*np.sin(theta)/total_cross_section

    # We want to perform rejection sampling, so we want to find a gaussian in theta which completely covers the PDF
    # PDF has max values at phi = pi/2 so covering at that phi ensures covering at all phi

    # Finding max of KN_PDF to centre gaussian cover
    if relative_energy==1.0: # Precaculated for the common case
        mean = 0.806
        rejection_scale = 0.45 
    else:
        # Want to find peak of distribution to within 0.15 radians
        check_values_at = np.linspace(0,np.pi/2,20)
        max_index = np.argmax(KN_PDF(check_values_at,np.pi/2))
        mean = check_values_at[max_index] 
        # The scaling constant should be A such that A/sqrt(2pisigma^2) > KN_PDF(mean), set A/sqrt(2pisigma^2)=1.2*KN_PDF(mean)
        # Because mean is at peak.
        # sqrt(2pisigma^2)=1
        rejection_scale = 1.2*KN_PDF(mean,np.pi/2)

    sigma = np.sqrt(np.pi/2)
    
    # Create the gaussian PDF
    gaussian = sc.stats.norm(loc = mean, scale = sigma)

    # Perform rejection sampling
    # Sample using the gauusian cover, then 
    # Uniformly pick a p between the scaled gaussian distribution at that theta,phi and 0
    # Then if that p is less than the Klein-Nishima distribution at that theta,phi, then accept the theta,phi, else reject and try again.
    while True:
        theta_try = sigma*np.random.randn()+mean
        phi_try = np.random.uniform(0,2*np.pi)
        p = np.random.uniform(0,rejection_scale*gaussian.pdf(theta_try)/(2*np.pi)) # Max value is cover_PDF(theta)*cover_PDF(phi)
        if p<=KN_PDF(theta_try, phi_try):
            return (theta_try, phi_try)

# All photons are instances of the Gamma class.
class Gamma:
    def __init__(self, source_pos:list, source_direction:list, source_pol:list, creation_time:float):
        # Source coordinates and direction, angles in radians
        self.source_pos = source_pos
        self.source_direction = source_direction # Direction of propagation when created
        self.source_polarisation = source_pol
        self.creation_time = creation_time

        # Coordinates of the previous point of interaction/creation, and the current polarisation and direction.
        self.last_pos = source_pos
        self.direction = source_direction
        self.polarisation = source_pol

        # Other properties
        self.energy = 511 # KeV
        self.compton_happened = False
        self.compton_pos:list # Position where compton happened
        self.detection_happened = False
        self.detect_pos:list # Position where detection happened
        self.exists = True # Becomes False if photon is absorbed. Preferred over deleting the photon to avoid possible bugs.
     
    def frequency(self): # Hz
        return self.energy/hbar
    
    def wavelength(self): # m
        return c/self.frequency()

    def momentum(self): # eV/c
        return self.energy/c
    
    def setEnergy(self,energy):
        self.energy = energy

    def compton(self, compton_radius:float, compton_shape = 0, compton_width = compton_width):
        '''
        Controls the changes in a photon's energy, direction, polarisation when it undergoes scattering.
        Takes in parameters of the scattering surface, checks for intersection and returns None if no intersection, else returns intersection point.
        Radius is ring/sphere radius, width is ring width, shape decides ring/sphere/potentially others.
        '''
        # Energy relative to mass energy of electrons
        epsilon = self.energy/(511)

        self.compton_pos = intersection(start_loc = self.last_pos, direction = self.direction, radius = compton_radius, width = compton_width, shape=compton_shape) 
        if not self.compton_pos:
            return None # If intersection() returns None, the photon did not intersect with the compton ring and hence did not scatter.
        self.last_pos = self.compton_pos
        self.compton_happened = True

        '''
        I'm not sure the below implementation for the rotation of polarisation is correct.
        Currently I'm interpeting the 1968 paper "Polarization of Compton-Scattered Gamma Rays"
        to be saying that the polarisation rotates the same way as the direction vector.
        That does make sense, since what is actually rotating is the electric field's direction of oscillation
        and in the same way the magnetic field's direction of oscillation, and these together determine
        the direction of motion, so it is reasonable to expect all three to have the same rotation.
        The possible flaw in this argument is what if the magnetic field rotates some other way.
        To answer that, I need to better understand the scattering process, and whether QED has answers to my questions.
        '''
        # Sample the KN distribution to find angles by which the photon scatters.
        # delta_theta is the scattering angle with respect to direction of motion.
        # delta_phi is the scattering angle with respect to the polarisation vector in the plane normal to original direction of motion.
        delta_theta, delta_phi = Sampling_Klein_Nishima(epsilon)
        original_direction = copy.deepcopy(self.direction)
        original_polarisation = copy.deepcopy(self.polarisation)

        '''
        The rotations are essentially as if original_direction is the z-axis, original_polarisation is the x-axis,
        And we have a vector at spherical angles delta_theta, delta_phi in this coordinate system.
        We want to move the z-axis to this vector, and the x-axis should undergo the same rotation.
        '''
        # Rotate by the polar angle in the x-z plane
        self.direction = vec.rotate_about(self.direction,
                                          vec.cross_product(original_direction,
                                                            original_polarisation),
                                          delta_theta)
        self.polarisation = vec.rotate_about(self.polarisation,
                                            vec.cross_product(original_direction,
                                                            original_polarisation),
                                            delta_theta)

        # Rotate by the azimuthal angle in the old x-y plane
        self.direction = vec.rotate_about(self.direction,
                                          original_direction,
                                          delta_phi)
        self.polarisation = vec.rotate_about(self.polarisation,
                                             original_direction,
                                             delta_phi)

        # Energy loss due to compton scattering.
        self.energy = self.energy/(1+epsilon*(1-np.cos(delta_theta)))

        return self.compton_pos

    def detect(self,detector_radius, detector_shape=0, detector_width = detector_width):
        '''
        Takes in parameters of the detection surface, checks for intersection and returns None if no intersection, else returns intersection point.
        '''
        if not self.exists:
            return None

        # Does the line start_loc + direction*t intersect the detector surface?
        
        self.detect_pos = intersection(start_loc = self.last_pos, direction = self.direction,
                                       radius = detector_radius, width = detector_width,
                                       shape = detector_shape) # Default settings choose detector parameters
        if self.detect_pos:
            self.detection_happened = True
            self.exists = False # Detection intrinsically involves absorption.
        return self.detect_pos
    
    def distance_travelled(self):
        '''
        Calculates total distance travelled by a photon in its lifetime until the last interaction we are aware of.
        Used to calculate time of detection.
        Can be updated to take in a time and compare with creation time of photon to give more complete answer,
        But since this is only used to calculate time of detection right now, that's not necessary right now.
        '''
        distance = 0

        if self.compton_happened:
            distance += vec.magnitude(vec.add(self.compton_pos,(-1)*self.source_pos))
            distance += vec.magnitude(vec.add(self.last_pos,(-1)*self.compton_pos)) # Mostly a useless line
        else:
            distance += vec.magnitude(vec.add(self.last_pos,(-1)*self.source_pos)) # Again, mostly useless
        
        if self.detection_happened:
            distance += vec.magnitude(vec.add(self.detect_pos,(-1)*self.last_pos))

        return distance

    def absorb(self, absorbtion_radius:float, absorber_shape=0, absorber_width = absorber_width):
        '''
        Takes in parameters of the absorbing surface, checks for intersection and returns None if no intersection, else returns intersection point.
        The absorption position it calculates is only accurate if the absorber radius is decided probabilistically from throughout the absorbing region.
        But since there is no practical need for that yet, the code treats it as if all absorption of the region happens at its boundary.
        '''
        absorbtion_pos = intersection(start_loc = self.last_pos, direction = self.direction, radius = absorbtion_radius, width = absorber_width, shape=absorber_shape)
        if absorbtion_pos:
            self.last_pos = absorbtion_pos
            self.exists = False
        return absorbtion_pos


# Create n random photons
def random_photons(n,time_window=0,start_time=0, source_location = [0,0,0]):
    '''
    Create n pairs of photons with random directions and polarisations at a given location.
    Can add generation in random location from a given region later.
    Used to generate all the photons for a time window, so distributed their creation times uniformly throughout that window.
    '''
    photons = np.empty([2*n], dtype = Gamma)

    # Generate n random pairs of angles for the direction unit vector.
    directions_theta = np.pi*np.random.random(n)
    directions_phi = 2*np.pi*np.random.random(n)
    # Once given the direction, the polarisation is randomly generated by taking any vector normal to the direction
    # And rotating it by a random angle, which is polarisation_phi.
    polarisation_phi = 2*np.pi*np.random.random(n)

    times = time_window*np.random.random(n)+start_time

    for i in range(n):
        direction = vec.get_cartesian(rho=1,theta=directions_theta[i], phi=directions_phi[i])

        # Make perpendicular vector for polarisation
        polarisation = vec.cross_product(direction,hatz) # Cross product with z unit vector
        if vec.are_parallel(direction,hatz) or vec.are_antiparallel(direction,hatz): # Except if the direction is the z axis, then
            polarisation = vec.cross_product(direction,hatx) # Cross product with x unit vector
        
        # Rotate by random angle
        polarisation = vec.rotate_about(polarisation,direction, polarisation_phi[i])

        photons[2*i] = Gamma(source_location,
                             direction,
                             polarisation,
                             times[i])
        
        # Each event generates two photons back to back and with perpendicular polarisations
        photons[2*i+1] = Gamma(source_location,
                               -direction,
                               vec.rotate_about(polarisation,direction,np.pi/2),
                               times[i])

    return photons

def simulation(time_unit, time_window, total_time, source_activity_slice, compton_probability = compton_probability, return_times = False, source_location = [0,0,0], attenuation = 0, silent=False):
    '''
    Performs a simulation, given time parameters, the activity per time_unit, compton and absorption parameters and point source location.
    Attenuation is absorption coefficient in mm^(-1).
    return_times decides if the time taken per photon generation, absorption, compton scattering and detection is returned to the user.
    These are total over all windows.
    List mode data returned in the detected_photons list.
    '''
    if not silent:
        print(f"Simulation started, will take approximately {round(total_time*source_activity_slice*5.7173221111/10**5/0.8,2)} seconds.")

    n_windows = int(total_time//time_window) # If total time is not an integral multiple, only simulate all the windows fitting inside it.
    source_activity_window = source_activity_slice*time_window/time_unit
    detected_photons = []
    times = [0,0,0,0]
    
    # loop
    for window in range(n_windows):
        start = time.time()

        # Source activity - number of degenerations this window - must be an integer.
        # So if source_activity_window = 10.3, source_activity = 10 30% of the time and 11 70% of the time.
        source_activity = int(source_activity_window)
        if np.random.random(1)<(source_activity_window - source_activity): # The difference is the probability of one more event this window
            source_activity+=1
        
        photons = random_photons(source_activity, time_window,window*time_window, source_location)
        
        photon_gen_done = time.time()
        times[0] += photon_gen_done - start

        # Attenuation
        for i in range(2*source_activity):
            photon = photons[i]
            if np.random.random(1)<(1-np.exp(-attenuation*absorber_radius)):
                photon.absorb(absorber_radius)


        attenuation_done = time.time()
        times[1] += attenuation_done - photon_gen_done

        # Compton scattering
        # Because the distributions are exactly symmetrical, we can simply solve for only the first photon compton scattering 
        # with the probability of either of them scattering.
        for i in range(source_activity):

            scattering_photon = photons[2*i]
            total_compton_probability = compton_probability*(2-compton_probability) # Probability of either photon compton scattering

            # In case the first photon has been absorbed, choose the second and make probab normal
            if not scattering_photon.exists:
                scattering_photon = photons[2*i+1]
                total_compton_probability = compton_probability
                # If second was also absorbed, no compton scattering
                if not scattering_photon.exists:
                    break
            # If second was absorbed but first wasn't, fix probability even then
            if not photons[2*i+1].exists:
                total_compton_probability = compton_probability
            
            if np.random.random(1)<total_compton_probability:
                scattering_photon.compton(compton_radius)

        compton_done = time.time()
        times[2] += compton_done - attenuation_done

        # Detection
        for i in range(2*source_activity):
            
            photon = photons[i]

            detection_position = photon.detect(detector_radius,detector_width=detector_width)
            
            if detection_position: # If not detected, this is None
                # There are 360 detectors, so detector number is simply the greatest integer function of phi in degrees + 1
                detector_number = int(180/np.pi*vec.get_phi(detection_position))+1
                
                # Energy blurring
                det_energy = np.random.normal(photon.energy, np.sqrt(photon.energy/(8*np.log(2)))) 
                
                # Check if energy is within detection range
                if (detector_energy_cutoff_low < det_energy < detector_energy_cutoff_high):
                    # Time blurring
                    det_time = np.random.normal(photon.creation_time + photon.distance_travelled()/c,
                                                timing_resolution) # timing resolution is FWHM of time blurring.

                    # Store the detected photons sorted by time
                    # Format of detected photons is time, detector number, detected energy, polarisation vector, pair identity
                    # (Window number, pair number) uniquely identifies each true pair

                    num_already_detected = len(detected_photons)
                    if num_already_detected==0:
                        # If no photons detected yet, simply append this one to the list.
                        detected_photons.append([det_time, detector_number,det_energy,photon.polarisation,(window,i//2)])
                    # else go through the list of detected photons in reverse, checking where this one lies temporally.
                    for j in range(1,num_already_detected+1):
                        # If the jth photon from the end was detected after this one, go back till you find the one detected before this.
                        if det_time < detected_photons[-j][0]: 
                            continue
                        # If the 1st one from the end was detected before this one, append to the end.
                        elif j==1:
                            detected_photons.append([det_time, detector_number,det_energy,photon.polarisation, (window,i//2)])
                            break
                        # If any other jth photon detected before this, insert appropriately.
                        else:
                            detected_photons.insert(-j+1,[det_time, detector_number,det_energy,photon.polarisation, (window,i//2)])
                            break
        # Store all created photons for future comparison? So far the plan is no.

        detection_done = time.time()
        times[3] += detection_done - compton_done
    
    if return_times:
        return detected_photons, times
    return detected_photons

def pairing(list_mode_detections, pairing_window = 10, polarisation_correlation = False, polarisation_min_angle_deg = 80, compton_filter = True):
    '''
    Pairs events submitted in list mode form, which is the format [time, detector number, energy, polarisation vector, pair identity].
    The list mode data must be sorted by time.
    pairing_window is size of the coincidence pairing window in ns.
    compton_filter decides whether low energy photons are filtered out before pairing.
    polarisation_correlation decides whether pairing is done based on polarisations. When True, no multiplets, only pairs grouped.
    From a multiplet, pairs with cos(theta) closest to 0 are paired.
    polarisation_min_angle_deg is the min angle between the polarisations. For reference, cos(80) = 0.173.
    Returns pairs or multiplets of grouped events. 
    '''
    event_pairs = [] # Pairs (or multiplets) of events

    i = 0
    num_detections = len(list_mode_detections)
    while i < num_detections - 1:
        # Start new window
        window_events = [list_mode_detections[i]]

        # Get all events within pairing window
        i+=1 # So i is always at the index after the last one included in window_events

        #  Check if the ith event's detection time is within the window
        while list_mode_detections[i][0] < window_events[0][0] + pairing_window: 
            window_events.append(list_mode_detections[i]) # If yes, append to window events
            i+=1 # And move to the next event
            if i == num_detections: # If the list mode data is exhausted, stop.
                break
        
        # Filter the objects in the window for scattering
        if compton_filter:
            # Check if the energy of the event is within sqrt(E) of E, E = 511 KeV
            window_events = list(filter((lambda event: abs(event[2] - 511) < np.sqrt(511)),
                                        window_events))
        
        # Do polarisation filtering
        if polarisation_correlation:
            num_events = len(window_events)
            polarisation_pairs = []
            # Calculate the cosines of the angles between polarisation vectors for each pair of events
            for j in range(num_events):
                for k in range(j+1,num_events):
                    polarisation_cosine = vec.dot_product(window_events[j][3],
                                                        window_events[k][3])
                    # Filter by some minimum angle separation - hence a maximum in cos(angle)
                    if polarisation_cosine < np.cos(np.radians(polarisation_min_angle_deg)):
                        # Refer to events by their indices in window_events
                        # And so store (j,k) tuples in pol pairs list
                        polarisation_pairs.append(((j,k),polarisation_cosine)) 
            # Sort by polarisation_cosines in ascending
            polarisation_pairs = sorted(polarisation_pairs,key = (lambda x: x[1]))
            
            while polarisation_pairs: # We go through the list until it is empty
                chosen_pair = polarisation_pairs[0] # The first ((i,j),cosine) tuple left in the list
                # Append the pair [event_i, event_j] to event pairs
                event_pairs.append([window_events[chosen_pair[0][0]],
                                    window_events[chosen_pair[0][1]]])
                # Remove any pairs which also have events indexed by i or j
                for pair in polarisation_pairs:
                    if (chosen_pair[0][0] in pair[0]) or (chosen_pair[0][1] in pair[0]):
                        polarisation_pairs.remove(pair)
                # Thus we also lose the chosen pair, and the loop will terminate only when the list is empty
        else: # If polarisation correlation not used
            # If there are multiple events in the window, consider the multiplet, else discard
            if len(window_events)>1:
                event_pairs.append(window_events)
        
        # And then the loop starts another window with no overlap with this one.
    return event_pairs

def num_true_pairings(event_multiplets, strict = False, random_pair = False):
    '''
    Checks how many of the groupings in a set of multiplets are correct.
    Strict mode discards all multiplets as false pairings.
    random_pair picks a random pair from the multiplet and checks that, and if its correct the multiplet is counted as a true pairing, else not.
    Without either of them, if any of the pairs in the multiplet are a true pair count it as true.
    If both settings are on, strict takes priority. 
    '''
    true_pairings = 0
    if strict:
        for multiplet in event_multiplets:
            # The 5th item in an detected event list is the pair identity tuple,
            # if that is the same then the events originate from the same disintegration.
            if len(multiplet)==2 and multiplet[0][4] == multiplet[1][4]: 
                true_pairings+=1
    elif random_pair:
        for multiplet in event_multiplets:
            len_multiplet = len(multiplet)
            while True:
                i = random.randrange(len_multiplet)
                j = random.randrange(len_multiplet)
                if i!=j:
                    break
            if multiplet[i][4]==multiplet[j][4]:
                    true_pairings+=1
    else:
        for multiplet in event_multiplets:
            len_multiplet = len(multiplet)
            # Iterate over all pairs from the multiplet
            for i,j in generate_pairs_double_iterate(len_multiplet):
                if multiplet[i][4]==multiplet[j][4]:
                    true_pairings+=1
                    break
    return true_pairings

def generate_pairs_double_iterate(length):
    '''
    Given a number of things, returns the indices for all possible pairs that could be chosen from the things.
    '''
    pairs_of_indices = []
    for i in range(length):
        for j in range(i+1,length):
            pairs_of_indices.append((i,j))
    return pairs_of_indices

def scattering_detection_vectors_plots(number_of_photons, source_location = [0,0,0]):
    '''
    Makes a 3D plot of the direction vectors of a number of photons 
    at their creation, at their compton scattering and at their detection.
    Considers compton probability = 1, spherical scatterer and detector surfaces, no attenuation.
    Draws the scatterer and detector surfaces, which are small (scatterer has 40 mm radius) for better visualisation.
    Can set source location.
    '''
    photons = random_photons(number_of_photons, source_location=source_location)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # n_photons = len(photons)

    all_vectors = [[],[],[]]
    origins = [[],[],[]]
    small_radius = 40
    detector_radius = 1.5*small_radius
    for i in range(len(photons)):

        photon = photons[i]
        photon.compton(small_radius, compton_shape = 1)

        for i in range(3):
            all_vectors[i].append(photon.source_direction[i])
            origins[i].append(source_location[i])
            all_vectors[i].append(photon.direction[i])
            origins[i].append(photon.compton_pos[i])

        photon.detect(detector_radius, detector_shape = 1)

        for i in range(3):
            all_vectors[i].append(photon.direction[i])
            origins[i].append(photon.detect_pos[i])


    ax.quiver(*origins, *all_vectors,
            length=small_radius/5)

    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]

    colours = ['red','orange']
    radii = [small_radius, detector_radius]
    for i in range(len(radii)):
        radius = radii[i]
        x = radius*np.cos(u)*np.sin(v)
        y = radius*np.sin(u)*np.sin(v)
        z = radius*np.cos(v)
        ax.plot_wireframe(x, y, z, color = colours[i], alpha = 0.3)

        # ax.set_xlim3d(-8,8)
        # ax.set_ylim3d(-8,8)
        # ax.set_zlim3d(-8,8)
    plt.axis(False)
    plt.show()
    return plt.figure

if __name__=="__main__":
    # Simulation
    # Let source at (rho,theta) = (70 mm, 30º):
    source_rho = 10 # mm
    source_theta = 30 # degrees
    total_time = 10**6
    source_activity_slice = 0.8

    detected_photons, times = simulation(time_unit, 100, total_time, source_activity_slice, compton_probability=0.25,return_times=True,
                                          source_location=[source_rho*np.cos(np.radians(source_theta)),
                                                          source_rho*np.sin(np.radians(source_theta)),
                                                          0],
                                        attenuation = 0.005) # With this attenuation we lose 30% of the photons
    # At attenuation 0.0096 per mm (water), we lose 60% of the photons
    print("Time taken:")
    print("Photon Generation"," Attenuation","      Compton Scattering","  Detection",sep=", ")
    print(*times,sep=", ")
    energy_paired_events = pairing(detected_photons, time_window)


    print("Detections:",len(detected_photons))


    energy_paired_events = pairing(detected_photons, time_window)
    print("Energy-paired multiplets:",len(energy_paired_events),"of which true pairs",num_true_pairings(energy_paired_events),", but being strict",num_true_pairings(energy_paired_events,strict=True))

    
    polarisation_paired_events = pairing(detected_photons, time_window, polarisation_correlation=True, polarisation_min_angle_deg=80)
    print("Polarisation pairs:",len(polarisation_paired_events),"of which true pairs",num_true_pairings(polarisation_paired_events))


    if make_plots:
        scattering_detection_vectors_plots(20, source_location=[source_rho*np.cos(np.radians(source_theta)),
                                                          source_rho*np.sin(np.radians(source_theta)),
                                                          0])