"""
This file contains functions that assist in calculating spatial properties. For
example, they assist in spatial orientation and alignment of the structures.
"""

# Numerical
import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.optimize import least_squares

# Biological
from Bio import PDB

# Internal
from . import load



def align_structures(rel_pps, def_pps, al_chain_ids, common_res=None,
                     rel_dict=None, def_dict=None, target="atoms", deleted_atom=[]):
    """
    Perform spatial alignment of deformed polypeptides to reference polypep-
    tides for the given chains.
    """
    # Load relaxed/deformed chains to align
    rel_al_ch = load.choose_chains(rel_pps, al_chain_ids)
    def_al_ch = load.choose_chains(def_pps, al_chain_ids)

    if target == "atoms":
        # Load relaxed/deformed common atoms to align
        rel_al_atm, _ = load.atoms(rel_al_ch, common_res, rel_dict, deleted_atom)
        def_al_atm, _ = load.atoms(def_al_ch, common_res, def_dict, deleted_atom)

        # Load all deformed atoms
        def_all_atom, _ = load.atoms(def_pps, common_res, def_dict, deleted_atom)

    elif target == "CA":
        # Load relaxed/deformed common atoms to align
        rel_al_atm, _ = load.atoms_CA(rel_al_ch, common_res, rel_dict)
        def_al_atm, _ = load.atoms_CA(def_al_ch, common_res, def_dict)

        # Load all deformed atoms
        def_all_atom, _ = load.atoms_CA(def_pps, common_res, def_dict)

        
    # Load alignment tool, align selected atoms, and apply to all
    super_imposer = PDB.Superimposer()
    super_imposer.set_atoms(rel_al_atm, def_al_atm)
    super_imposer.apply(def_all_atom)
    
    return # need to specify which atoms, otherwise can not perform alignment!

def align_structures_from_index(rel_pps, def_pps, al_idx, 
                                common_res=None, rel_dict=None, 
                                def_dict=None, deleted_atom=[]):
    """
    Perform spatial alignment of deformed polypeptides to reference polypep-
    tides for the given indexes.
    """
    # Load all relaxed/deformed atoms
    rel_all_atom, _ = load.atoms_CA(rel_pps, common_res, rel_dict, deleted_atom)
    def_all_atom, _ = load.atoms_CA(def_pps, common_res, def_dict, deleted_atom)
    
    # Load relaxed/deformed common atoms to align
    rel_al_atm = [rel_all_atom[i] for i in al_idx]
    def_al_atm = [def_all_atom[i] for i in al_idx]
    
    # Load alignment tool, align selected atoms, and apply to all
    super_imposer = PDB.Superimposer()
    super_imposer.set_atoms(rel_al_atm, def_al_atm)
    super_imposer.apply(def_all_atom)
    
    return # need to specify which atoms, otherwise can not perform alignment!

def cylinder_axis(xyz):

    """
    Fit the provided coordinates to a cylinder of variable center, orientation
    and radius (5 parameters). Returns orientation axis of the best fit
    cylinder, as well as the whole fit results.
    """   
    # Initialize fit parameters
    xc, yc = np.mean(xyz[:,0]), np.mean(xyz[:,1]) # cylinder center
    theta, phi = 0., .1 # angles about x and y axis, phi>0 biases Z>0
    r = np.std(xyz) # radius

    # Least square fit using cylinder distance
    result = least_squares(cylinder_error,
                           [xc, yc, theta, phi, r],
                           args=(xyz[:,0], xyz[:,1], xyz[:,2]),
                           max_nfev=10000)
    [xc, yc, theta, phi, r] = result.x

    # Calculate new axis
    r = R.from_euler('xy', [theta, phi])
    axis = r.apply([0, 0, 1])

    return axis, result


def cylinder_error(p, x, y, z):
    """
    Error from the data to a cylindrical sheet centerd in xc, yc = p[0], p[1]
    of radius r = p[4], and with x/y orientation theta = p[2] and phi = p[3].
    """
    # Create xyz array
    xyz = np.transpose(np.array([x, y, z]))
    
    # Position vectors to cylinder center
    dxyz = xyz - np.array([p[0], p[1], 0.]) 
    
    # Cylinder axis
    r = R.from_euler('xy', [p[2], p[3]])
    cyl_axis = r.apply([0, 0, 1])
    
    # Distance to axis
    dxyz_projected = dxyz * cyl_axis
    axis_distance = np.sqrt(np.sum(dxyz**2, 1) - np.sum(dxyz_projected, 1)**2)
    fit_error = np.sum((axis_distance - p[4])**2)
    
    return fit_error


def center_and_align(pps, cyl_chain_ids, Nrep=2):
    """
    This function aligns the given polypeptides to the  z axis by using the
    chains in cyl_chain_ids as a cylinder that will be centerd and aligned. We
    repeat the mean subtraction, as after a finite rotation the center gets
    displaced.
    """
    
    # Load coordinates of cylinder chains
    cyl_pps = load.choose_chains(pps, cyl_chain_ids)
    cyl_xyz, _ = load.coordinates(cyl_pps)
    
    # Calculate displacement and axis
    displacement = np.mean(cyl_xyz, 0)
    axis, _ = cylinder_axis(cyl_xyz - displacement)
    rigid_body(pps, displacement, axis)
    
    # Center again
    #cyl_pps = load.choose_chains(pps, cyl_chain_ids)
    cyl_xyz, _ = load.coordinates(cyl_pps)
    displacement = np.mean(cyl_xyz, 0)
    rigid_body(pps, displacement)

    return axis


def rigid_body(pps, displacement, old_axis=[0,0,1], new_axis=[0,0,1]):
    """
    Move polypeptides by displacement and rotate them so that the
    input axis matches the Z axis.
    """
    # Define rotation matrix
    rotation = PDB.rotmat(PDB.Vector(new_axis),
                          PDB.Vector(old_axis))
    
    # Perform rotation of residues
    for chains in pps:
        for pp in chains:
            for res in pp:
                if load.test_residue(res):
                    res.transform(rotation,
                                  -displacement)
        
    return

def rigid_body_general(chains, displacement, old_axis=[0,0,1], new_axis=[0,0,1]):
    """
    Move chains by displacement and rotate them so that the
    input axis matches the Z axis.
    """
    # Define rotation matrix
    rotation = PDB.rotmat(PDB.Vector(new_axis),
                          PDB.Vector(old_axis))
    
    # Perform rotation of residues
    for chain in chains:
        for res in chain:
            for at in res:
                at.transform(rotation, -displacement)

    return


def noise(xyz, sigma):
    """
    asd
    """
    
    # Create noise data structure
    xyz_n = np.zeros_like(xyz)

    # Add normal noise
    xyz_n[:, :] = xyz[:, :] + np.random.normal(0, sigma, (len(xyz), 3))

    return xyz_n


def generate_rod(r0=20., z0=180., theta0=2*np.pi, ds = 4.):
    """
    Generate a cylindrical rod with the provided dimensions (Radius height, in
    A), and the provided spacing among atoms (ds)
    """

    # Then, assign coordiante values
    rtz = []
    Nr0, Nz0 = int(r0 / ds), int(z0 / ds)

    for z in np.linspace(0, z0, Nz0) :
        for r in np.linspace(0, r0, Nr0):
            Ntheta0 = int(2 * np.pi * r / ds)
            for theta in np.linspace(0, theta0, Ntheta0):
                rtz.append([r, theta, z])

    # Transform to cartesian
    rtz = np.array(rtz)
    xyz = coordinates_to_polar(rtz, inv=True)

    return xyz

def generate_sphere(r0=20., layers = 5, density = 5):
    """
    Generate a sphere provided dimensions (radius, in A), and the provided 
    density
    """

    rtp = []
    rs = np.arange(1, layers+1)
    n = rs**2*density
    
    for i in range(layers):
        for j in range(n[i]):
            r = r0*i/layers
            t = np.random.uniform(0, np.pi)
            p = np.random.uniform(0, 2*np.pi)
            rtp.append([r, t, p])
    
    rtp = np.array(rtp)
    xyz = coordinates_to_spherical(rtp, True)
    
    return xyz


def deform(xyz, method='twist', d=0.01):
    """
    Take as input a set of coordinates, a method, and a deformation parameter;
    and select the corresponding deformation function
    """
    cases = {'twist': twist,
             'spin': spin,
             'radial': radial_extension,
             'shear': shear}

    return cases[method](xyz, d)


def twist(xyz, d=0.01):
    """
    Apply a twist transformation to the given (euclidean) coordinates
    """
    # Transform reference coordinates to cylindrical
    rtz = coordinates_to_polar(xyz)
    rtz_d = []

    # Apply twist transformation
    for rtz_a in rtz:
        r_d = rtz_a[0]
        t_d = rtz_a[1] + d * rtz_a[2]
        z_d = rtz_a[2]
        rtz_d.append([r_d, t_d, z_d])

    # Transform deformed coordinates to euclidean
    rtz_d = np.array(rtz_d)
    xyz_d = coordinates_to_polar(rtz_d, inv=True)

    return xyz_d


def radial_extension(xyz, d=.5):
    """
    Apply a uniform radial extension to the given (euclidean) coordinates
    """
    # Transform reference coordinates to cylindrical
    rtz = coordinates_to_polar(xyz)
    rtz_d = []

    # Apply radaial extensiomn transformation
    for rtz_a in rtz:
        r_d = rtz_a[0] + d * rtz_a[0]
        t_d = rtz_a[1] 
        z_d = rtz_a[2]
        rtz_d.append([r_d, t_d, z_d])

    # Transform deformed coordinates to euclidean
    rtz_d = np.array(rtz_d)
    xyz_d = coordinates_to_polar(rtz_d, inv=True)

    return xyz_d


def spin(xyz, d):
    """
    Deform angles in the radial direction
    """
    # Transform reference coordinates cylindrical
    rtz = coordinates_to_polar(xyz)
    rtz_d = []

    # Apply spin
    for rtz_a in rtz:
        r_d = rtz_a[0]
        t_d = rtz_a[1] + d * rtz_a[0]
        z_d = rtz_a[2]
        rtz_d.append([r_d, t_d, z_d])

    # Transform deformed coordinates to euclidean
    rtz_d = np.array(rtz_d)
    xyz_d = coordinates_to_polar(rtz_d, inv=True)

    return xyz_d


def shear(xyz, d=0.1):
    """
    Apply a shear transformation to the given (euclidean) coordinates
    """
    xyz_d = []
    
    # Apply shear transformation
    for xyz_a in xyz:
        x_d = xyz_a[0] + d * xyz_a[2]
        y_d = xyz_a[1] + d * xyz_a[2]
        z_d = xyz_a[2]
        xyz_d.append([x_d, y_d, z_d])

    xyz_d = np.array(xyz_d)
    
    return xyz_d


def coordinates_to_polar(coordinates, inv=False):
    """
    Transform the given coordinates from cartesian to polar, or from polar to
    cartesian (inverse) if inv=True
    """
    if not inv:
        xyz = coordinates
        r = np.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2)
        t = np.arctan2(xyz[:, 1], xyz[:, 0])
        z = xyz[:, 2]
        rtz =  np.transpose(np.array([r, t, z]))
        return rtz
    else:
        rtz = coordinates
        x = rtz[:, 0] * np.cos(rtz[:, 1])
        y = rtz[:, 0] * np.sin(rtz[:, 1])
        z = rtz[:, 2]
        xyz = np.transpose(np.array([x, y, z]))
        return xyz


def tensor_to_polar(xyz, T, inv=False):
    """
    Transform the given tensor to polar coordinates using the given
    (cartesian) coordinates.
    """
    # Obtain polar angle
    rtz = coordinates_to_polar(xyz)
    theta = rtz[:, 1]
    
    # Define rotation matrix and its transpose
    R = np.array([np.array([ np.cos(theta),
                  np.sin(theta),
                  np.zeros_like(theta)]),
                  np.array([-np.sin(theta),
                  np.cos(theta),
                  np.zeros_like(theta)]),
                  np.array([np.zeros_like(theta),
                  np.zeros_like(theta),
                  np.ones_like(theta)])])
    R = np.transpose(R, (2,0,1))
    Rt = np.transpose(R, (0,2,1))
    
    # Perform dot product
    R_x_T = np.einsum('abc,acf->abf', R, T)
    T_cyl = np.einsum('abc,acf->abf', R_x_T, Rt)
    
    return T_cyl


def vector_to_polar(xyz, v, inv=False):
    """
    Transform the given vector to polar coordinates using the given (cartesian)
    coordinates.
    """
    # Obtain polar angle
    rtz = coordinates_to_polar(xyz)
    theta = rtz[:, 1]
    
    # Define rotation matrix and its transpose
    R = np.array([np.array([ np.cos(theta),
                  np.sin(theta),
                  np.zeros_like(theta)]),
                  np.array([-np.sin(theta),
                  np.cos(theta),
                  np.zeros_like(theta)]),
                  np.array([np.zeros_like(theta),
                  np.zeros_like(theta),
                  np.ones_like(theta)])])
    R = np.transpose(R, (2, 0, 1))

    # Perform dot product
    R_x_v = np.einsum('aij,aj->ai', R, v)

    return R_x_v


def coordinates_to_spherical(coordinates, inv=False):
    """
    Transform the given coordinates from cartesian to spherical, or from spherical to
    cartesian (inverse) if inv=True
    """
    if not inv:
        xyz = coordinates
        xy2 = xyz[:,0]**2 + xyz[:,1]**2
        
        r = np.sqrt(xy2 + xyz[:,2]**2)
        t = np.arctan2(np.sqrt(xy2), xyz[:,2])
        p = np.arctan2(xyz[:,1], xyz[:,0])
        rtp =  np.transpose(np.array([r, t, p]))
        return rtp
    else:
        rtp = coordinates
        x = rtp[:, 0] * np.sin(rtp[:, 1]) * np.cos(rtp[:, 2])
        y = rtp[:, 0] * np.sin(rtp[:, 1]) * np.sin(rtp[:, 2])
        z = rtp[:, 0] * np.cos(rtp[:, 1])
        xyz = np.transpose(np.array([x, y, z]))
        return xyz