"""
Galaxy Calculation Utilities - Astropy Ready
============================================

Comprehensive galaxy physics calculations with full Astropy integration.
Includes morphology, dynamics, photometry, and evolution calculations.
"""

import logging
from typing import Dict, Union
import numpy as np

# Astropy imports - direct imports following standards
import astropy.units as u
import astropy.constants as const
import astropy.coordinates as coord
from astropy.cosmology import Planck18, WMAP9, FlatLambdaCDM
from astropy.coordinates import SkyCoord

# Import from our astropy integration
from .astropy_integration import (
    AstronomicalConstants, CoordinateSystem, UnitConversions, 
    CosmologyCalculations, quick_distance
)

logger = logging.getLogger(__name__)

# Use Astropy constants
c = const.c
H0_default = 70 * u.km / u.s / u.Mpc  # Default Hubble constant
M_sun = const.M_sun
L_sun = const.L_sun
pc = const.pc


def calculate_galaxy_luminosity_distance(redshift: float, 
                                       cosmology: str = 'Planck18') -> u.Quantity:
    """
    Calculate luminosity distance to galaxy using Astropy cosmology.

    Args:
        redshift: Galaxy redshift
        cosmology: Cosmology model ('Planck18', 'WMAP9', 'FlatLambdaCDM')

    Returns:
        Luminosity distance with units (Mpc)
    """
    cosmo_calc = CosmologyCalculations(cosmology)
    d_L = cosmo_calc.luminosity_distance(redshift)
    return d_L * u.Mpc


def calculate_galaxy_angular_size(physical_size: u.Quantity, 
                                distance: u.Quantity) -> u.Quantity:
    """
    Calculate angular size of galaxy using Astropy units.

    Args:
        physical_size: Physical size with units (kpc, pc, etc.)
        distance: Distance with units (Mpc, kpc, etc.)

    Returns:
        Angular size with units (arcsec)
    """
    angular_size = UnitConversions.physical_to_angular_size(
        physical_size.to(u.kpc).value,
        distance.to(u.Mpc).value
    )
    return angular_size * u.arcsec


def calculate_galaxy_distance_modulus(distance: u.Quantity) -> u.Quantity:
    """
    Calculate distance modulus for galaxy using Astropy units.

    Args:
        distance: Distance with units

    Returns:
        Distance modulus with units (mag)
    """
    distance_pc = distance.to(u.pc)
    mu = 5 * np.log10(distance_pc.value) - 5
    return mu * u.mag


def calculate_galaxy_mass_from_luminosity(luminosity: u.Quantity, 
                                        mass_to_light_ratio: u.Quantity = 2.0*u.M_sun/u.L_sun
                                        ) -> u.Quantity:
    """
    Calculate galaxy mass from luminosity using Astropy units.

    Args:
        luminosity: Luminosity with units (L_sun, erg/s, etc.)
        mass_to_light_ratio: Mass-to-light ratio with units

    Returns:
        Galaxy mass with units (M_sun)
    """
    luminosity_solar = luminosity.to(u.L_sun)
    ml_ratio_solar = mass_to_light_ratio.to(u.M_sun / u.L_sun)
    
    mass = luminosity_solar * ml_ratio_solar
    return mass.to(u.M_sun)


def calculate_galaxy_rotation_curve(radius: u.Quantity, total_mass: u.Quantity,
                                   disk_scale_length: u.Quantity = 3.0*u.kpc
                                   ) -> u.Quantity:
    """
    Calculate galaxy rotation curve using Astropy units.

    Args:
        radius: Radii with units
        total_mass: Total galaxy mass with units
        disk_scale_length: Disk scale length with units

    Returns:
        Rotation velocities with units (km/s)
    """
    # Convert to standard units
    r = radius.to(u.kpc)
    M_tot = total_mass.to(u.M_sun)
    R_d = disk_scale_length.to(u.kpc)
    
    # Component mass fractions
    disk_fraction = 0.3
    halo_fraction = 0.7
    
    # Disk component (exponential profile)
    M_disk = M_tot * disk_fraction
    x = r / R_d
    M_disk_enclosed = M_disk * (1 - np.exp(-x) * (1 + x))
    
    # Halo component (NFW profile approximation)
    M_halo = M_tot * halo_fraction
    R_s = 20 * u.kpc  # Typical scale radius
    M_halo_enclosed = M_halo * (r / (r + R_s))
    
    # Total enclosed mass
    M_enclosed = M_disk_enclosed + M_halo_enclosed
    
    # Circular velocity
    v_circ = np.sqrt(const.G * M_enclosed / r)
    return v_circ.to(u.km / u.s)


def calculate_galaxy_density_profile(radius: u.Quantity, 
                                    central_density: u.Quantity,
                                    scale_length: u.Quantity,
                                    profile_type: str = "exponential"
                                    ) -> u.Quantity:
    """
    Calculate galaxy surface density profile using Astropy units.

    Args:
        radius: Radii with units
        central_density: Central surface density with units (M_sun/pc²)
        scale_length: Scale length with units
        profile_type: "exponential", "sersic", or "deVaucouleurs"

    Returns:
        Surface density with units (M_sun/pc²)
    """
    r = radius.to(u.kpc)
    Sigma_0 = central_density.to(u.M_sun / u.pc**2)
    R_s = scale_length.to(u.kpc)
    
    if profile_type == "exponential":
        # Exponential disk profile
        density = Sigma_0 * np.exp(-r / R_s)
        
    elif profile_type == "sersic":
        # Sersic profile (n=2 default)
        n = 2.0  # Sersic index
        b_n = 1.9992 * n - 0.3271  # Approximation for b_n
        density = Sigma_0 * np.exp(-b_n * ((r / R_s)**(1/n) - 1))
        
    elif profile_type == "deVaucouleurs":
        # de Vaucouleurs profile (Sersic n=4)
        b_4 = 7.6692
        density = Sigma_0 * np.exp(-b_4 * ((r / R_s)**0.25 - 1))
        
    else:
        raise ValueError(f"Unknown profile type: {profile_type}")
    
    return density


def calculate_star_formation_rate(galaxy_mass: u.Quantity, 
                                 galaxy_type: str = "Sb", 
                                 redshift: float = 0.0) -> u.Quantity:
    """
    Calculate star formation rate using empirical relations with Astropy units.

    Args:
        galaxy_mass: Galaxy stellar mass with units
        galaxy_type: Galaxy morphological type
        redshift: Galaxy redshift

    Returns:
        Star formation rate with units (M_sun/year)
    """
    mass_solar = galaxy_mass.to(u.M_sun).value
    log_mass = np.log10(mass_solar)
    
    # Main sequence relation (log-linear)
    if log_mass < 9.5:
        log_sfr_base = 0.8 * log_mass - 7.2
    else:
        log_sfr_base = 0.3 * log_mass - 2.5
    
    # Morphology correction
    morph_correction = {
        "E0": -1.5, "E3": -1.2, "E7": -1.0, "S0": -0.8,
        "Sa": -0.3, "Sb": 0.0, "Sc": 0.3,
        "SBa": -0.2, "SBb": 0.1, "SBc": 0.4,
        "Irr": 0.6,
    }
    
    correction = morph_correction.get(galaxy_type, 0.0)
    
    # Redshift evolution (more star formation at higher z)
    z_evolution = 0.5 * redshift
    
    log_sfr = log_sfr_base + correction + z_evolution
    sfr = 10**log_sfr
    
    return sfr * u.M_sun / u.year


def calculate_galaxy_color_index(stellar_mass: u.Quantity, 
                                star_formation_rate: u.Quantity,
                                metallicity: float = 0.02) -> u.Quantity:
    """
    Calculate galaxy color index (g-r) from physical properties using Astropy units.

    Args:
        stellar_mass: Stellar mass with units
        star_formation_rate: Star formation rate with units
        metallicity: Metallicity (fraction of solar)

    Returns:
        g-r color index with units (mag)
    """
    M_star = stellar_mass.to(u.M_sun).value
    SFR = star_formation_rate.to(u.M_sun / u.year).value
    
    # Specific star formation rate
    log_ssfr = np.log10(SFR / M_star)
    
    # Mass dependence (more massive galaxies are redder)
    log_mass = np.log10(M_star)
    mass_term = 0.15 * (log_mass - 10.0)
    
    # Star formation dependence (more star formation = bluer)
    sfr_term = -0.25 * (log_ssfr + 10.0)
    
    # Metallicity dependence (higher metallicity = redder)
    metal_term = 5.0 * (metallicity - 0.02)
    
    # Base color for typical galaxy
    base_color = 0.65
    
    color_gr = base_color + mass_term + sfr_term + metal_term
    color_gr = np.clip(color_gr, 0.2, 1.5)
    
    return color_gr * u.mag


def calculate_galaxy_metallicity(stellar_mass: u.Quantity, 
                                star_formation_rate: u.Quantity) -> float:
    """
    Calculate galaxy metallicity using mass-metallicity relation with Astropy units.

    Args:
        stellar_mass: Stellar mass with units
        star_formation_rate: Star formation rate with units

    Returns:
        Metallicity as fraction of solar (Z/Z_sun)
    """
    log_mass = np.log10(stellar_mass.to(u.M_sun).value)
    log_sfr = np.log10(star_formation_rate.to(u.M_sun / u.year).value)
    
    # Fundamental metallicity relation (Mannucci et al. 2010)
    mu = log_mass - 0.32 * log_sfr
    
    if mu < 9.5:
        metallicity_solar = 0.1 * 10**(0.6 * mu - 5.4)
    else:
        metallicity_solar = 0.6 * 10**(0.3 * mu - 2.85)
    
    return metallicity_solar


def calculate_galaxy_age_from_color(color_gr: u.Quantity, 
                                   metallicity: float = 1.0) -> u.Quantity:
    """
    Estimate galaxy age from color index using Astropy units.

    Args:
        color_gr: g-r color index with units
        metallicity: Metallicity relative to solar

    Returns:
        Age with units (Gyr)
    """
    color_value = color_gr.to(u.mag).value
    
    # Empirical relation from stellar population models
    base_age = 2.0 + 8.0 * (color_value - 0.3)
    
    # Metallicity correction (higher metallicity = younger for given color)
    metal_correction = -1.0 * np.log10(metallicity)
    
    age_gyr = base_age + metal_correction
    age_gyr = np.clip(age_gyr, 0.1, 13.8)  # Age of universe
    
    return age_gyr * u.Gyr


def calculate_galaxy_size_from_mass(stellar_mass: u.Quantity, 
                                   galaxy_type: str = "Sb") -> u.Quantity:
    """
    Calculate galaxy effective radius from stellar mass using Astropy units.

    Args:
        stellar_mass: Stellar mass with units
        galaxy_type: Galaxy morphological type

    Returns:
        Effective radius with units (kpc)
    """
    log_mass = np.log10(stellar_mass.to(u.M_sun).value)
    
    # Mass-size relation (different for early/late types)
    if galaxy_type in ["E0", "E3", "E7", "S0"]:
        # Early-type galaxies
        log_re = 0.56 * log_mass - 4.54
    else:
        # Late-type galaxies
        log_re = 0.22 * log_mass - 1.43
    
    effective_radius = 10**log_re
    return effective_radius * u.kpc


def calculate_galaxy_surface_brightness(luminosity: u.Quantity, 
                                       effective_radius: u.Quantity,
                                       distance: u.Quantity,
                                       band: str = "r") -> u.Quantity:
    """
    Calculate galaxy surface brightness using Astropy units.

    Args:
        luminosity: Luminosity with units
        effective_radius: Effective radius with units
        distance: Distance to galaxy with units
        band: Photometric band

    Returns:
        Surface brightness with units (mag/arcsec²)
    """
    # Convert to angular effective radius
    re_angular = calculate_galaxy_angular_size(effective_radius, distance)
    
    # Effective area in square arcseconds
    area_arcsec2 = np.pi * re_angular.to(u.arcsec).value**2
    
    # Convert luminosity to apparent magnitude
    luminosity_solar = luminosity.to(u.L_sun).value
    
    # Distance modulus
    mu = calculate_galaxy_distance_modulus(distance).value
    
    # Absolute magnitude (assuming solar absolute magnitude)
    M_abs = -2.5 * np.log10(luminosity_solar) + 4.83  # Solar abs mag in r-band
    
    # Apparent magnitude
    m_app = M_abs + mu
    
    # Surface brightness
    mu_surface = m_app + 2.5 * np.log10(area_arcsec2)
    
    return mu_surface * u.mag / u.arcsec**2


def calculate_bulge_to_disk_ratio(galaxy_type: str, stellar_mass: u.Quantity) -> float:
    """
    Calculate bulge-to-disk mass ratio using Astropy units.

    Args:
        galaxy_type: Galaxy morphological type
        stellar_mass: Total stellar mass with units

    Returns:
        Bulge-to-disk ratio (dimensionless)
    """
    # Base ratios by type
    base_ratios = {
        "E0": 100.0, "E3": 100.0, "E7": 100.0,  # Pure bulge
        "S0": 1.5, "Sa": 0.8, "Sb": 0.3, "Sc": 0.1,
        "SBa": 0.6, "SBb": 0.25, "SBc": 0.08,
        "Irr": 0.02,
    }
    
    base_ratio = base_ratios.get(galaxy_type, 0.3)
    
    # Mass dependence (more massive galaxies have larger bulges)
    log_mass = np.log10(stellar_mass.to(u.M_sun).value)
    mass_correction = 10**(0.2 * (log_mass - 10.5))
    
    return base_ratio * mass_correction


def calculate_galaxy_virial_mass(stellar_mass: u.Quantity, 
                                galaxy_type: str = "Sb") -> u.Quantity:
    """
    Calculate galaxy virial (total) mass from stellar mass using Astropy units.

    Args:
        stellar_mass: Stellar mass with units
        galaxy_type: Galaxy morphological type

    Returns:
        Virial mass with units (M_sun)
    """
    # Stellar-to-halo mass relation (abundance matching)
    log_M_star = np.log10(stellar_mass.to(u.M_sun).value)
    
    # Different relations for different types
    if galaxy_type in ["E0", "E3", "E7"]:
        # Early-type galaxies (higher stellar fractions)
        log_M_halo = 1.15 * log_M_star + 1.4
    elif galaxy_type == "S0":
        log_M_halo = 1.1 * log_M_star + 1.2
    else:
        # Late-type galaxies (lower stellar fractions)
        log_M_halo = 1.2 * log_M_star + 0.8
    
    M_virial = 10**log_M_halo
    return M_virial * u.M_sun


def calculate_galaxy_tidal_radius(galaxy_mass: u.Quantity, 
                                 neighbor_mass: u.Quantity,
                                 separation: u.Quantity) -> u.Quantity:
    """
    Calculate tidal radius for galaxy interactions using Astropy units.

    Args:
        galaxy_mass: Galaxy mass with units
        neighbor_mass: Neighboring galaxy mass with units  
        separation: Separation distance with units

    Returns:
        Tidal radius with units (kpc)
    """
    # Tidal radius calculation
    mass_ratio = galaxy_mass / neighbor_mass
    r_tidal = separation * (mass_ratio / 3)**(1/3)
    
    return r_tidal.to(u.kpc)


def get_galaxy_properties(galaxy_type: str, stellar_mass: u.Quantity, 
                         redshift: float = 0.0, distance: u.Quantity = None
                         ) -> Dict[str, Union[u.Quantity, float]]:
    """
    Calculate comprehensive galaxy properties using Astropy units.

    Args:
        galaxy_type: Galaxy morphological type
        stellar_mass: Stellar mass with units
        redshift: Galaxy redshift
        distance: Distance with units (calculated from redshift if None)

    Returns:
        Dict containing all galaxy properties with proper units
    """
    # Calculate distance if not provided
    if distance is None and redshift > 0:
        distance = calculate_galaxy_luminosity_distance(redshift)
    elif distance is None:
        distance = 10.0 * u.Mpc  # Default distance
    
    # Calculate derived properties
    sfr = calculate_star_formation_rate(stellar_mass, galaxy_type, redshift)
    metallicity = calculate_galaxy_metallicity(stellar_mass, sfr)
    color_gr = calculate_galaxy_color_index(stellar_mass, sfr, metallicity)
    age = calculate_galaxy_age_from_color(color_gr, metallicity)
    size = calculate_galaxy_size_from_mass(stellar_mass, galaxy_type)
    virial_mass = calculate_galaxy_virial_mass(stellar_mass, galaxy_type)
    
    # Angular size
    angular_size = calculate_galaxy_angular_size(size, distance)
    
    # Surface brightness
    luminosity = stellar_mass / (2.0 * u.M_sun / u.L_sun)  # Rough M/L ratio
    surface_brightness = calculate_galaxy_surface_brightness(luminosity, size, distance)
    
    # Bulge-to-disk ratio
    bulge_disk_ratio = calculate_bulge_to_disk_ratio(galaxy_type, stellar_mass)
    
    # Distance modulus
    distance_modulus = calculate_galaxy_distance_modulus(distance)
    
    return {
        "stellar_mass": stellar_mass.to(u.M_sun),
        "virial_mass": virial_mass,
        "star_formation_rate": sfr,
        "metallicity": metallicity,
        "color_gr": color_gr,
        "age": age,
        "effective_radius": size,
        "distance": distance,
        "redshift": redshift,
        "angular_size": angular_size,
        "surface_brightness": surface_brightness,
        "bulge_disk_ratio": bulge_disk_ratio,
        "luminosity": luminosity,
        "distance_modulus": distance_modulus,
    }


def get_galaxy_morphology_params(galaxy_type: str) -> Dict[str, Union[float, bool]]:
    """
    Get morphological parameters for galaxy type.

    Args:
        galaxy_type: Galaxy morphological type

    Returns:
        Dict containing morphological parameters
    """
    from .astronomical_data import GALAXY_TYPES
    
    if galaxy_type not in GALAXY_TYPES:
        logger.warning(f"Unknown galaxy type: {galaxy_type}, using Sb")
        galaxy_type = "Sb"
    
    params = GALAXY_TYPES[galaxy_type].copy()
    
    # Add derived parameters
    params["is_early_type"] = galaxy_type in ["E0", "E3", "E7", "S0"]
    params["has_spiral_arms"] = galaxy_type in ["Sa", "Sb", "Sc", "SBa", "SBb", "SBc"]
    params["has_bar"] = galaxy_type.startswith("SB")
    params["is_irregular"] = galaxy_type == "Irr"
    
    return params


# Legacy functions for backward compatibility
def calculate_galaxy_luminosity_distance_legacy(redshift: float, H0: float = 70.0, 
                                               omega_m: float = 0.3, 
                                               omega_lambda: float = 0.7) -> float:
    """Legacy function - use calculate_galaxy_luminosity_distance for new code"""
    import warnings
    warnings.warn("Legacy function - use calculate_galaxy_luminosity_distance with Astropy units", 
                  DeprecationWarning, stacklevel=2)
    
    result = calculate_galaxy_luminosity_distance(redshift, 'Planck18')
    return result.to(u.Mpc).value


def calculate_galaxy_angular_size_legacy(physical_size_kpc: float, 
                                        distance_mpc: float) -> float:
    """Legacy function - use calculate_galaxy_angular_size for new code"""
    import warnings
    warnings.warn("Legacy function - use calculate_galaxy_angular_size with Astropy units", 
                  DeprecationWarning, stacklevel=2)
    
    result = calculate_galaxy_angular_size(physical_size_kpc * u.kpc, distance_mpc * u.Mpc)
    return result.to(u.arcsec).value
