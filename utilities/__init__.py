"""
Astronomical Utilities - Astropy Ready
=====================================

Comprehensive utility functions for astronomical calculations and data processing.
Fully integrated with Astropy for units, constants, coordinates, and cosmology.
"""

# Astronomical data utilities
from .astronomical_data import (
    GALAXY_TYPES,
    HR_DIAGRAM_PARAMS,
    STELLAR_CLASSIFICATION,
    create_sample_stellar_data,
    get_galaxy_config,
    get_stellar_data,
    validate_hr_diagram_data,
    create_sample_hr_diagram_data,
    get_stellar_color_from_temperature,
    calculate_absolute_magnitude,
    calculate_distance_modulus,
    distance_modulus_to_parsecs,
)

# Astropy integration utilities
from .astropy_integration import (
    AstronomicalConstants,
    CoordinateSystem,
    AstronomicalTime,
    WorldCoordinateSystem,
    AstropyIO,
    CosmologyCalculations,
    UnitConversions,
    TableOperations,
    quick_skycoord,
    quick_time,
    quick_distance,
    quick_wcs_from_fits,
)

# FITS utilities (AstroPhot integration)
from .fits_utilities import (
    validate_fits_file,
    get_fits_info,
    suggest_visualization_params,
    prepare_fits_for_blender,
    generate_astronomical_coordinates,
    estimate_object_type,
    create_fits_metadata,
    calculate_surface_brightness,
    create_false_color_rgb,
    get_filter_color,
    FITS_EXTENSIONS,
    COMMON_FILTERS,
)

# Galaxy calculation utilities (Astropy ready)
from .galaxy_utilities import (
    calculate_galaxy_luminosity_distance,
    calculate_galaxy_angular_size,
    calculate_galaxy_distance_modulus,
    calculate_galaxy_mass_from_luminosity,
    calculate_galaxy_rotation_curve,
    calculate_galaxy_density_profile,
    calculate_star_formation_rate,
    calculate_galaxy_color_index,
    calculate_galaxy_metallicity,
    calculate_galaxy_age_from_color,
    calculate_galaxy_size_from_mass,
    calculate_galaxy_surface_brightness,
    calculate_bulge_to_disk_ratio,
    calculate_galaxy_virial_mass,
    calculate_galaxy_tidal_radius,
    get_galaxy_properties,
    get_galaxy_morphology_params,
    # Legacy functions (with deprecation warnings)
    calculate_galaxy_luminosity_distance_legacy,
    calculate_galaxy_angular_size_legacy,
)

# Physics calculation utilities (Astropy ready)
from .physics_utilities import (
    calculate_gravitational_force,
    calculate_escape_velocity,
    calculate_orbital_period,
    calculate_orbital_velocity,
    calculate_orbital_elements,
    calculate_orbital_position,
    calculate_hill_sphere,
    calculate_roche_limit,
    calculate_tidal_force,
    calculate_hohmann_transfer,
    calculate_gravitational_potential,
    calculate_virial_velocity,
    calculate_dark_matter_profile,
    propagate_orbit_astropy,
    create_astropy_orbit,
    # Legacy functions (with deprecation warnings)
    propagate_orbit,
    create_poliastro_orbit,
)

# Nebula calculation utilities
from .nebula_utilities import (
    calculate_nebula_angular_size,
    calculate_nebula_color_index,
    calculate_nebula_density,
    calculate_nebula_distance_modulus,
    calculate_nebula_expansion_velocity,
    calculate_nebula_luminosity,
    calculate_nebula_mass,
    calculate_nebula_metallicity,
    calculate_nebula_size,
    calculate_nebula_temperature,
    get_nebula_properties,
    calculate_nebula_ionization_parameter,
    calculate_emission_line_ratios,
)

# Post-processing utilities
from .post_processing_utilities import (
    calculate_bloom_intensity,
    calculate_chromatic_aberration,
    calculate_color_temperature_to_rgb,
    calculate_contrast_curve,
    calculate_depth_of_field,
    calculate_exposure_value,
    calculate_film_grain_intensity,
    calculate_lens_distortion,
    calculate_motion_blur_angle,
    calculate_saturation_multiplier,
    calculate_vignette_intensity,
    get_cinematic_preset,
    get_dramatic_preset,
    get_dreamy_preset,
    get_scientific_preset,
    apply_astrophotography_processing,
    calculate_signal_to_noise_ratio,
)

# Export all available functionality
__all__ = [
    # Astronomical data
    "GALAXY_TYPES",
    "HR_DIAGRAM_PARAMS", 
    "STELLAR_CLASSIFICATION",
    "create_sample_stellar_data",
    "get_galaxy_config",
    "get_stellar_data",
    "validate_hr_diagram_data",
    "create_sample_hr_diagram_data",
    "get_stellar_color_from_temperature",
    "calculate_absolute_magnitude",
    "calculate_distance_modulus",
    "distance_modulus_to_parsecs",
    
    # Astropy integration
    "AstronomicalConstants",
    "CoordinateSystem",
    "AstronomicalTime", 
    "WorldCoordinateSystem",
    "AstropyIO",
    "CosmologyCalculations",
    "UnitConversions",
    "TableOperations",
    "quick_skycoord",
    "quick_time",
    "quick_distance",
    "quick_wcs_from_fits",
    
    # FITS utilities
    "validate_fits_file",
    "get_fits_info",
    "suggest_visualization_params",
    "prepare_fits_for_blender",
    "generate_astronomical_coordinates",
    "estimate_object_type",
    "create_fits_metadata",
    "calculate_surface_brightness",
    "create_false_color_rgb",
    "get_filter_color",
    "FITS_EXTENSIONS",
    "COMMON_FILTERS",
    
    # Galaxy utilities (Astropy ready)
    "calculate_galaxy_luminosity_distance",
    "calculate_galaxy_angular_size",
    "calculate_galaxy_distance_modulus",
    "calculate_galaxy_mass_from_luminosity",
    "calculate_galaxy_rotation_curve",
    "calculate_galaxy_density_profile",
    "calculate_star_formation_rate",
    "calculate_galaxy_color_index",
    "calculate_galaxy_metallicity",
    "calculate_galaxy_age_from_color",
    "calculate_galaxy_size_from_mass",
    "calculate_galaxy_surface_brightness",
    "calculate_bulge_to_disk_ratio",
    "calculate_galaxy_virial_mass",
    "calculate_galaxy_tidal_radius",
    "get_galaxy_properties",
    "get_galaxy_morphology_params",
    
    # Physics utilities (Astropy ready)
    "calculate_gravitational_force",
    "calculate_escape_velocity",
    "calculate_orbital_period",
    "calculate_orbital_velocity",
    "calculate_orbital_elements",
    "calculate_orbital_position",
    "calculate_hill_sphere",
    "calculate_roche_limit",
    "calculate_tidal_force",
    "calculate_hohmann_transfer",
    "calculate_gravitational_potential",
    "calculate_virial_velocity",
    "calculate_dark_matter_profile",
    "propagate_orbit_astropy",
    "create_astropy_orbit",
    
    # Nebula utilities
    "calculate_nebula_angular_size",
    "calculate_nebula_color_index",
    "calculate_nebula_density",
    "calculate_nebula_distance_modulus",
    "calculate_nebula_expansion_velocity",
    "calculate_nebula_luminosity",
    "calculate_nebula_mass",
    "calculate_nebula_metallicity",
    "calculate_nebula_size",
    "calculate_nebula_temperature",
    "get_nebula_properties",
    "calculate_nebula_ionization_parameter",
    "calculate_emission_line_ratios",
    
    # Post-processing utilities
    "calculate_bloom_intensity",
    "calculate_chromatic_aberration",
    "calculate_color_temperature_to_rgb",
    "calculate_contrast_curve",
    "calculate_depth_of_field",
    "calculate_exposure_value",
    "calculate_film_grain_intensity",
    "calculate_lens_distortion",
    "calculate_motion_blur_angle",
    "calculate_saturation_multiplier",
    "calculate_vignette_intensity",
    "get_cinematic_preset",
    "get_dramatic_preset",
    "get_dreamy_preset",
    "get_scientific_preset",
    "apply_astrophotography_processing",
    "calculate_signal_to_noise_ratio",
]

# Astropy integration status
ASTROPY_AVAILABLE = True
ASTROPY_VERSION = None

try:
    import astropy
    ASTROPY_VERSION = astropy.__version__
    print(f"✅ Astropy {ASTROPY_VERSION} integration active")
except ImportError:
    ASTROPY_AVAILABLE = False
    print("⚠️ Astropy not available - some functions will have limited functionality")

# Quick access to common constants with Astropy units
try:
    from astropy import units as u, constants as const
    
    # Commonly used constants
    CONSTANTS = {
        'c': const.c,           # Speed of light
        'G': const.G,           # Gravitational constant  
        'h': const.h,           # Planck constant
        'k_B': const.k_B,       # Boltzmann constant
        'M_sun': const.M_sun,   # Solar mass
        'R_sun': const.R_sun,   # Solar radius
        'L_sun': const.L_sun,   # Solar luminosity
        'au': const.au,         # Astronomical unit
        'pc': const.pc,         # Parsec
    }
    
    # Commonly used units
    UNITS = {
        # Length
        'm': u.m, 'km': u.km, 'pc': u.pc, 'kpc': u.kpc, 'Mpc': u.Mpc,
        'au': u.au, 'ly': u.lyr, 'Rsun': u.R_sun,
        # Mass  
        'kg': u.kg, 'Msun': u.M_sun, 'Mearth': u.M_earth,
        # Time
        's': u.s, 'min': u.min, 'hour': u.hour, 'day': u.day, 
        'year': u.year, 'Myr': u.Myr, 'Gyr': u.Gyr,
        # Energy/Luminosity
        'J': u.J, 'erg': u.erg, 'eV': u.eV, 'Lsun': u.L_sun,
        # Angular
        'deg': u.deg, 'rad': u.rad, 'arcsec': u.arcsec, 'arcmin': u.arcmin,
        # Spectral
        'Hz': u.Hz, 'nm': u.nm, 'um': u.um, 'AA': u.AA,
        # Magnitudes
        'mag': u.mag,
    }
    
except ImportError:
    CONSTANTS = {}
    UNITS = {}
    print("⚠️ Astropy constants and units not available")

def get_available_utilities():
    """Get list of available utility modules and their Astropy status"""
    status = {
        'astronomical_data': {'available': True, 'astropy_ready': False},
        'astropy_integration': {'available': ASTROPY_AVAILABLE, 'astropy_ready': True},
        'fits_utilities': {'available': True, 'astropy_ready': True},
        'galaxy_utilities': {'available': True, 'astropy_ready': ASTROPY_AVAILABLE},
        'physics_utilities': {'available': True, 'astropy_ready': ASTROPY_AVAILABLE},
        'nebula_utilities': {'available': True, 'astropy_ready': False},  # TODO: Update
        'post_processing_utilities': {'available': True, 'astropy_ready': False},  # TODO: Update
    }
    
    return status

def check_astropy_dependencies():
    """Check which Astropy submodules are available"""
    dependencies = {}
    
    if not ASTROPY_AVAILABLE:
        return dependencies
    
    astropy_modules = [
        'astropy.units', 'astropy.constants', 'astropy.coordinates',
        'astropy.time', 'astropy.io.fits', 'astropy.io.ascii',
        'astropy.wcs', 'astropy.table', 'astropy.cosmology'
    ]
    
    for module in astropy_modules:
        try:
            __import__(module)
            dependencies[module] = True
        except ImportError:
            dependencies[module] = False
    
    return dependencies

def create_astropy_example_usage():
    """Create examples showing Astropy integration"""
    if not ASTROPY_AVAILABLE:
        return "Astropy not available - install with: pip install astropy"
    
    examples = """
# Astropy Integration Examples
# ===========================

from utilities import *
import astropy.units as u

# 1. Galaxy calculations with units
stellar_mass = 1e11 * u.M_sun
distance = 50 * u.Mpc
galaxy_props = get_galaxy_properties('Sb', stellar_mass)

# 2. Coordinate transformations
coord = quick_skycoord(150.0, -30.0)  # RA, Dec in degrees
galactic = coord.galactic

# 3. Physics with proper units
orbit_data = create_astropy_orbit(
    semi_major_axis=1.0 * u.au,
    eccentricity=0.1,
    central_mass=1.0 * u.M_sun
)

# 4. Cosmological calculations
redshift = 0.1
distance = quick_distance(redshift)

# 5. FITS file handling
wcs = quick_wcs_from_fits('galaxy.fits')
"""
    return examples

# Print integration status on import
if __name__ != "__main__":
    status = get_available_utilities()
    astropy_ready_count = sum(1 for s in status.values() if s['astropy_ready'])
    total_count = len(status)
    
    print(f"📊 Astropy integration: {astropy_ready_count}/{total_count} modules ready")
    
    if ASTROPY_AVAILABLE:
        deps = check_astropy_dependencies()
        failed_deps = [k for k, v in deps.items() if not v]
        if failed_deps:
            print(f"⚠️ Missing Astropy modules: {', '.join(failed_deps)}")
        else:
            print("✅ All Astropy dependencies available")
    else:
        print("💡 Install Astropy for full functionality: pip install astropy")
