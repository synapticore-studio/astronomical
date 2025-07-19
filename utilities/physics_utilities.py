"""
Physics Calculation Utilities - Astropy Ready
=============================================

Comprehensive physics calculations for astronomical objects and systems.
Fully integrated with Astropy for units, constants, and coordinate systems.
"""

import numpy as np
from typing import Dict, Any, Tuple, Union
import logging

# Astropy imports - direct imports following standards
import astropy.units as u
import astropy.constants as const
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.time import Time

# Import from our astropy integration
from .astropy_integration import AstronomicalConstants, CoordinateSystem, UnitConversions

logger = logging.getLogger(__name__)

# Use Astropy constants
G = const.G
c = const.c
M_sun = const.M_sun
M_earth = const.M_earth
R_sun = const.R_sun
R_earth = const.R_earth
au = const.au
pc = const.pc
L_sun = const.L_sun


def calculate_gravitational_force(mass1: u.Quantity, mass2: u.Quantity, 
                                distance: u.Quantity) -> u.Quantity:
    """
    Calculate gravitational force between two masses using Astropy units.
    
    Args:
        mass1: Mass of first object with units (e.g., M_sun, kg)
        mass2: Mass of second object with units
        distance: Distance between objects with units (e.g., au, pc, m)
        
    Returns:
        Gravitational force with units (N)
    """
    force = G * mass1 * mass2 / (distance**2)
    return force.to(u.N)


def calculate_escape_velocity(mass: u.Quantity, radius: u.Quantity) -> u.Quantity:
    """
    Calculate escape velocity from an object using Astropy units.
    
    Args:
        mass: Object mass with units
        radius: Object radius with units
        
    Returns:
        Escape velocity with units (km/s)
    """
    v_esc = np.sqrt(2 * G * mass / radius)
    return v_esc.to(u.km / u.s)


def calculate_orbital_period(semi_major_axis: u.Quantity, 
                           central_mass: u.Quantity) -> u.Quantity:
    """
    Calculate orbital period using Kepler's third law with Astropy units.
    
    Args:
        semi_major_axis: Semi-major axis with units
        central_mass: Central mass with units
        
    Returns:
        Orbital period with units (years)
    """
    period = 2 * np.pi * np.sqrt(semi_major_axis**3 / (G * central_mass))
    return period.to(u.year)


def calculate_orbital_velocity(orbital_radius: u.Quantity, 
                             central_mass: u.Quantity) -> u.Quantity:
    """
    Calculate circular orbital velocity using Astropy units.
    
    Args:
        orbital_radius: Orbital radius with units
        central_mass: Central mass with units
        
    Returns:
        Orbital velocity with units (km/s)
    """
    v_orb = np.sqrt(G * central_mass / orbital_radius)
    return v_orb.to(u.km / u.s)


def calculate_orbital_elements(position: u.Quantity, velocity: u.Quantity,
                              central_mass: u.Quantity) -> Dict[str, u.Quantity]:
    """
    Calculate orbital elements from position and velocity vectors with Astropy units.
    
    Args:
        position: Position vector with units (3D)
        velocity: Velocity vector with units (3D)
        central_mass: Central mass with units
        
    Returns:
        Dict containing orbital elements with proper units
    """
    # Ensure inputs have correct dimensions
    if position.ndim == 1:
        position = position.reshape(1, -1)
    if velocity.ndim == 1:
        velocity = velocity.reshape(1, -1)
    
    r = np.linalg.norm(position, axis=-1)
    v = np.linalg.norm(velocity, axis=-1)
    
    # Specific orbital energy
    mu = G * central_mass
    energy = v**2 / 2 - mu / r
    
    # Semi-major axis
    a = -mu / (2 * energy)
    
    # Angular momentum vector
    h_vec = np.cross(position, velocity, axis=-1)
    h = np.linalg.norm(h_vec, axis=-1)
    
    # Eccentricity vector
    r_unit = position / r[..., np.newaxis]
    e_vec = np.cross(velocity, h_vec, axis=-1) / mu.value - r_unit
    e = np.linalg.norm(e_vec, axis=-1)
    
    # Inclination
    i = np.arccos(h_vec[..., 2] / h) * u.rad
    
    # Orbital period
    period = calculate_orbital_period(a, central_mass)
    
    return {
        "semi_major_axis": a.to(u.au),
        "eccentricity": e * u.dimensionless_unscaled,
        "inclination": i.to(u.deg),
        "period": period,
        "energy": energy.to(u.J / u.kg),
        "angular_momentum": h.to(u.m**2 / u.s)
    }


def calculate_orbital_position(semi_major_axis: u.Quantity, eccentricity: float,
                              true_anomaly: u.Quantity, inclination: u.Quantity = 0*u.deg,
                              longitude_ascending_node: u.Quantity = 0*u.deg,
                              argument_of_periapsis: u.Quantity = 0*u.deg) -> u.Quantity:
    """
    Calculate position from orbital elements using Astropy units.
    
    Args:
        semi_major_axis: Semi-major axis with units
        eccentricity: Orbital eccentricity (dimensionless)
        true_anomaly: True anomaly with angular units
        inclination: Inclination with angular units
        longitude_ascending_node: Longitude of ascending node with angular units
        argument_of_periapsis: Argument of periapsis with angular units
        
    Returns:
        Position vector with units matching semi_major_axis
    """
    # Convert angles to radians
    nu = true_anomaly.to(u.rad).value
    i = inclination.to(u.rad).value
    omega = longitude_ascending_node.to(u.rad).value
    w = argument_of_periapsis.to(u.rad).value
    
    # Distance from focus
    r = semi_major_axis * (1 - eccentricity**2) / (1 + eccentricity * np.cos(nu))
    
    # Position in orbital plane
    x_orb = r * np.cos(nu)
    y_orb = r * np.sin(nu)
    z_orb = 0 * r.unit
    
    # Rotation matrices
    cos_w, sin_w = np.cos(w), np.sin(w)
    cos_i, sin_i = np.cos(i), np.sin(i) 
    cos_omega, sin_omega = np.cos(omega), np.sin(omega)
    
    # Combined rotation matrix elements
    R11 = cos_w * cos_omega - sin_w * cos_i * sin_omega
    R12 = -sin_w * cos_omega - cos_w * cos_i * sin_omega
    R21 = cos_w * sin_omega + sin_w * cos_i * cos_omega
    R22 = -sin_w * sin_omega + cos_w * cos_i * cos_omega
    R31 = sin_w * sin_i
    R32 = cos_w * sin_i
    
    # Transform to inertial frame
    x = R11 * x_orb + R12 * y_orb
    y = R21 * x_orb + R22 * y_orb
    z = R31 * x_orb + R32 * y_orb
    
    return np.array([x, y, z]) * r.unit


def calculate_hill_sphere(primary_mass: u.Quantity, secondary_mass: u.Quantity,
                         orbital_distance: u.Quantity) -> u.Quantity:
    """
    Calculate Hill sphere radius for gravitational influence using Astropy units.
    
    Args:
        primary_mass: Primary mass with units
        secondary_mass: Secondary mass with units
        orbital_distance: Orbital separation with units
        
    Returns:
        Hill sphere radius with same units as orbital_distance
    """
    mass_ratio = secondary_mass / primary_mass
    hill_radius = orbital_distance * (mass_ratio / 3)**(1/3)
    return hill_radius


def calculate_roche_limit(primary_mass: u.Quantity, primary_radius: u.Quantity,
                         secondary_density: u.Quantity) -> u.Quantity:
    """
    Calculate Roche limit for tidal disruption using Astropy units.
    
    Args:
        primary_mass: Primary mass with units
        primary_radius: Primary radius with units
        secondary_density: Secondary density with units
        
    Returns:
        Roche limit distance with same units as primary_radius
    """
    # Calculate primary density
    primary_volume = (4/3) * np.pi * primary_radius**3
    primary_density = primary_mass / primary_volume
    
    # Rigid body Roche limit
    roche_limit = 2.44 * primary_radius * (primary_density / secondary_density)**(1/3)
    return roche_limit


def calculate_tidal_force(primary_mass: u.Quantity, distance: u.Quantity,
                         object_size: u.Quantity) -> u.Quantity:
    """
    Calculate tidal force on extended object using Astropy units.
    
    Args:
        primary_mass: Primary mass with units
        distance: Distance to primary with units
        object_size: Size of object experiencing tidal force with units
        
    Returns:
        Tidal acceleration with units (m/s²)
    """
    tidal_accel = 2 * G * primary_mass * object_size / (distance**3)
    return tidal_accel.to(u.m / u.s**2)


def calculate_hohmann_transfer(r1: u.Quantity, r2: u.Quantity, 
                              central_mass: u.Quantity) -> Dict[str, u.Quantity]:
    """
    Calculate Hohmann transfer orbit parameters using Astropy units.
    
    Args:
        r1: Initial orbital radius with units
        r2: Final orbital radius with units
        central_mass: Central mass with units
        
    Returns:
        Dict containing transfer orbit parameters with proper units
    """
    # Transfer orbit semi-major axis
    a_transfer = (r1 + r2) / 2
    
    # Velocities
    v1_circular = calculate_orbital_velocity(r1, central_mass)
    v2_circular = calculate_orbital_velocity(r2, central_mass)
    
    # Transfer velocities
    mu = G * central_mass
    v1_transfer = np.sqrt(mu * (2/r1 - 1/a_transfer))
    v2_transfer = np.sqrt(mu * (2/r2 - 1/a_transfer))
    
    # Delta-v requirements
    delta_v1 = abs(v1_transfer - v1_circular)
    delta_v2 = abs(v2_circular - v2_transfer)
    delta_v_total = delta_v1 + delta_v2
    
    # Transfer time
    transfer_period = calculate_orbital_period(a_transfer, central_mass)
    transfer_time = transfer_period / 2
    
    return {
        "transfer_semi_major_axis": a_transfer.to(u.au),
        "delta_v1": delta_v1.to(u.km/u.s),
        "delta_v2": delta_v2.to(u.km/u.s),
        "total_delta_v": delta_v_total.to(u.km/u.s),
        "transfer_time": transfer_time.to(u.year),
        "transfer_period": transfer_period.to(u.year)
    }


def calculate_gravitational_potential(mass: u.Quantity, radius: u.Quantity) -> u.Quantity:
    """
    Calculate gravitational potential at surface using Astropy units.
    
    Args:
        mass: Object mass with units
        radius: Object radius with units
        
    Returns:
        Gravitational potential with units (J/kg)
    """
    potential = -G * mass / radius
    return potential.to(u.J / u.kg)


def calculate_virial_velocity(mass: u.Quantity, radius: u.Quantity) -> u.Quantity:
    """
    Calculate virial velocity for gravitational system using Astropy units.
    
    Args:
        mass: System mass with units
        radius: System radius with units
        
    Returns:
        Virial velocity with units (km/s)
    """
    v_virial = np.sqrt(G * mass / radius)
    return v_virial.to(u.km / u.s)


def calculate_dark_matter_profile(radius: u.Quantity, mass_200: u.Quantity,
                                 concentration: float = 10.0,
                                 profile_type: str = "NFW") -> u.Quantity:
    """
    Calculate dark matter halo density profile using Astropy units.
    
    Args:
        radius: Radii with units (kpc, Mpc, etc.)
        mass_200: Virial mass with units (M_sun)
        concentration: Halo concentration parameter
        profile_type: "NFW", "Einasto", or "Burkert"
        
    Returns:
        Density profile with units (M_sun/pc³)
    """
    # Convert inputs to standard units
    r = radius.to(u.pc)
    M_200 = mass_200.to(u.M_sun)
    
    # Critical density (Planck 2018 cosmology)
    H0 = 67.4 * u.km / u.s / u.Mpc  # Hubble constant
    rho_crit = 3 * H0**2 / (8 * np.pi * G)
    rho_crit = rho_crit.to(u.M_sun / u.pc**3)
    
    # Assume Omega_m = 0.315 for matter density
    Omega_m = 0.315
    rho_m = Omega_m * rho_crit
    
    # Virial radius
    Delta_vir = 200  # Overdensity factor
    r_200 = (3 * M_200 / (4 * np.pi * Delta_vir * rho_m))**(1/3)
    
    # Scale radius
    r_s = r_200 / concentration
    
    if profile_type == "NFW":
        # Navarro-Frenk-White profile
        x = r / r_s
        
        # Characteristic density
        delta_c = (Delta_vir/3) * concentration**3 / (np.log(1 + concentration) - concentration/(1 + concentration))
        rho_s = delta_c * rho_m
        
        # NFW profile
        density = rho_s / (x * (1 + x)**2)
        
    elif profile_type == "Einasto":
        # Einasto profile
        alpha = 0.17  # Typical value
        r_e = r_s  # Use scale radius as effective radius
        
        # Normalization (approximate)
        from scipy.special import gamma
        rho_e = M_200 / (4 * np.pi * r_e**3 * alpha * gamma(3/alpha) * np.exp(2/alpha))
        
        # Einasto profile
        density = rho_e * np.exp(-2/alpha * ((r/r_e)**alpha - 1))
        
    elif profile_type == "Burkert":
        # Burkert profile (cored)
        r_b = r_s
        # Normalization factor for Burkert profile
        norm_factor = 4 * np.pi * r_b**3 * (1 + np.log(2))
        rho_b = M_200 / norm_factor
        
        # Burkert profile
        density = rho_b / ((1 + r/r_b) * (1 + (r/r_b)**2))
        
    else:
        raise ValueError(f"Unknown profile type: {profile_type}")
    
    return density.to(u.M_sun / u.pc**3)


def propagate_orbit_astropy(initial_position: u.Quantity, initial_velocity: u.Quantity,
                           central_mass: u.Quantity, time_array: u.Quantity,
                           method: str = "rk4") -> Tuple[u.Quantity, u.Quantity]:
    """
    Propagate orbital motion using numerical integration with Astropy units.
    
    Args:
        initial_position: Initial position vector with units (3D)
        initial_velocity: Initial velocity vector with units (3D)
        central_mass: Central mass with units
        time_array: Array of time points with units
        method: Integration method ("euler", "rk4")
        
    Returns:
        Tuple of (positions, velocities) arrays with proper units
    """
    # Convert to base units for integration
    pos_unit = initial_position.unit
    vel_unit = initial_velocity.unit
    time_unit = time_array.unit
    
    pos = initial_position.to(u.m).value
    vel = initial_velocity.to(u.m/u.s).value
    mass = central_mass.to(u.kg).value
    times = time_array.to(u.s).value
    
    n_steps = len(times)
    positions = np.zeros((n_steps, 3))
    velocities = np.zeros((n_steps, 3))
    
    # Initial conditions
    positions[0] = pos
    velocities[0] = vel
    
    # Gravitational parameter in SI units
    mu = G.to(u.m**3 / u.kg / u.s**2).value * mass
    
    for i in range(1, n_steps):
        dt = times[i] - times[i-1]
        
        if method == "euler":
            # Simple Euler integration
            r = np.linalg.norm(positions[i-1])
            acceleration = -mu * positions[i-1] / r**3
            
            velocities[i] = velocities[i-1] + acceleration * dt
            positions[i] = positions[i-1] + velocities[i-1] * dt
            
        elif method == "rk4":
            # Runge-Kutta 4th order
            def derivatives(pos_vec, vel_vec):
                r = np.linalg.norm(pos_vec)
                acc = -mu * pos_vec / r**3
                return vel_vec, acc
            
            # RK4 integration
            pos_curr = positions[i-1]
            vel_curr = velocities[i-1]
            
            k1_pos, k1_vel = derivatives(pos_curr, vel_curr)
            k2_pos, k2_vel = derivatives(pos_curr + k1_pos*dt/2, vel_curr + k1_vel*dt/2)
            k3_pos, k3_vel = derivatives(pos_curr + k2_pos*dt/2, vel_curr + k2_vel*dt/2)
            k4_pos, k4_vel = derivatives(pos_curr + k3_pos*dt, vel_curr + k3_vel*dt)
            
            positions[i] = pos_curr + dt/6 * (k1_pos + 2*k2_pos + 2*k3_pos + k4_pos)
            velocities[i] = vel_curr + dt/6 * (k1_vel + 2*k2_vel + 2*k3_vel + k4_vel)
    
    # Convert back to original units
    positions_with_units = positions * u.m
    velocities_with_units = velocities * u.m / u.s
    
    return positions_with_units.to(pos_unit), velocities_with_units.to(vel_unit)


def create_astropy_orbit(semi_major_axis: u.Quantity, eccentricity: float,
                        inclination: u.Quantity = 0*u.deg,
                        central_mass: u.Quantity = 1*u.M_sun) -> Dict[str, u.Quantity]:
    """
    Create orbit parameters with full Astropy unit support.
    
    Args:
        semi_major_axis: Semi-major axis with units
        eccentricity: Orbital eccentricity
        inclination: Inclination with angular units
        central_mass: Central mass with units
        
    Returns:
        Dict containing comprehensive orbit parameters with units
    """
    # Orbital period
    period = calculate_orbital_period(semi_major_axis, central_mass)
    
    # Mean motion
    mean_motion = 2 * np.pi * u.rad / period
    
    # Orbital velocities
    v_perihelion = calculate_orbital_velocity(semi_major_axis * (1 - eccentricity), central_mass)
    v_aphelion = calculate_orbital_velocity(semi_major_axis * (1 + eccentricity), central_mass)
    v_mean = calculate_orbital_velocity(semi_major_axis, central_mass)
    
    # Specific orbital energy
    specific_energy = -G * central_mass / (2 * semi_major_axis)
    
    # Angular momentum per unit mass
    h = np.sqrt(G * central_mass * semi_major_axis * (1 - eccentricity**2))
    
    return {
        "semi_major_axis": semi_major_axis.to(u.au),
        "eccentricity": eccentricity * u.dimensionless_unscaled,
        "inclination": inclination.to(u.deg),
        "central_mass": central_mass.to(u.M_sun),
        "period": period.to(u.year),
        "mean_motion": mean_motion.to(u.deg/u.day),
        "velocity_perihelion": v_perihelion.to(u.km/u.s),
        "velocity_aphelion": v_aphelion.to(u.km/u.s),
        "velocity_mean": v_mean.to(u.km/u.s),
        "specific_energy": specific_energy.to(u.J/u.kg),
        "angular_momentum": h.to(u.m**2/u.s),
        "perihelion_distance": (semi_major_axis * (1 - eccentricity)).to(u.au),
        "aphelion_distance": (semi_major_axis * (1 + eccentricity)).to(u.au)
    }


# Legacy functions for backward compatibility (with deprecation warnings)
def propagate_orbit(*args, **kwargs):
    """Legacy function - use propagate_orbit_astropy for new code"""
    import warnings
    warnings.warn("propagate_orbit is deprecated, use propagate_orbit_astropy", 
                  DeprecationWarning, stacklevel=2)
    # Convert to Astropy and call new function
    return propagate_orbit_astropy(*args, **kwargs)

def create_poliastro_orbit(*args, **kwargs):
    """Legacy function - use create_astropy_orbit for new code"""
    import warnings
    warnings.warn("create_poliastro_orbit is deprecated, use create_astropy_orbit", 
                  DeprecationWarning, stacklevel=2)
    return create_astropy_orbit(*args, **kwargs)
