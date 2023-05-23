import numpy as np
from scipy.spatial.transform import Rotation as R

def add(vec1,vec2):
    return [x+y for (x,y) in zip(vec1,vec2)]

def cross_product(vec1, vec2):
    return np.array([vec1[(i+1)%3]*vec2[(i+2)%3]-vec1[(i+2)%3]*vec2[(i+1)%3] for i in range(3)])

def dot_product(vec1, vec2):
    return sum([vec1[i]*vec2[i] for i in range(3)])

# Polar coordinates (Assuming cartesian input)
def get_r(vec):
    vec = np.array(vec)
    return np.sqrt(sum(vec**2))

def magnitude(vec):
    return get_r(vec)

def get_theta(vec):
    return np.arccos(vec[2]/get_r(vec))

def get_phi(vec):
    # Edge cases
    if vec[0]==0:
        if vec[1]==0:
            return 0
        else:
            return np.pi/2
    # Normal function
    return (np.arctan2(vec[0],vec[1])) % (2*np.pi)
    # except:
    #     # print(np)
    #     print(vec[0],vec[1])
    #     return (np.arctan2(vec[0],vec[1])) % (2*np.pi)

def get_cartesian(rho,theta,phi):
    return float(rho)*np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])

def rotate_about(vec, axis, angle):
    r = R.from_rotvec(angle*np.array(axis)/get_r(axis))
    return r.apply(vec)

def are_parallel(vec1,vec2):
    return (abs(dot_product(vec1,vec2) - magnitude(vec1)*magnitude(vec2)) < 1e-5)

def are_antiparallel(vec1,vec2):
    return (abs(dot_product(vec1,vec2) + magnitude(vec1)*magnitude(vec2)) < 1e-5)

