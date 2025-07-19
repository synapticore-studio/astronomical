"""
Astronomical Data Utilities
===========================

Comprehensive astronomical data definitions, sample data generators, and validation functions.
Self-contained implementation without external dependencies.
"""

import numpy as np
from typing import Dict, List, Any, Optional, Union

# Stellar Classification System (Morgan-Keenan)
STELLAR_CLASSIFICATION = {
    "O": {
        "temperature_range": (30000, 50000),
        "mass_range": (15, 90),
        "luminosity_range": (30000, 1000000),
        "color": (0.6, 0.7, 1.0),
        "lifetime_myr": 10,
        "main_sequence_fraction": 0.00003,
        "examples": ["Alnitak", "Mintaka", "Alnilam"]
    },
    "B": {
        "temperature_range": (10000, 30000),
        "mass_range": (2.1, 16),
        "luminosity_range": (25, 30000),
        "color": (0.7, 0.8, 1.0),
        "lifetime_myr": 400,
        "main_sequence_fraction": 0.0013,
        "examples": ["Rigel", "Spica", "Regulus"]
    },
    "A": {
        "temperature_range": (7500, 10000),
        "mass_range": (1.4, 2.1),
        "luminosity_range": (5, 25),
        "color": (0.9, 0.9, 1.0),
        "lifetime_myr": 2000,
        "main_sequence_fraction": 0.006,
        "examples": ["Sirius", "Vega", "Altair"]
    },
    "F": {
        "temperature_range": (6000, 7500),
        "mass_range": (1.04, 1.4),
        "luminosity_range": (1.5, 5),
        "color": (1.0, 1.0, 0.9),
        "lifetime_myr": 7000,
        "main_sequence_fraction": 0.03,
        "examples": ["Procyon", "Canopus", "Polaris"]
    },
    "G": {
        "temperature_range": (5200, 6000),
        "mass_range": (0.8, 1.04),
        "luminosity_range": (0.6, 1.5),
        "color": (1.0, 1.0, 0.8),
        "lifetime_myr": 12000,
        "main_sequence_fraction": 0.076,
        "examples": ["Sun", "Alpha Centauri A", "Capella"]
    },
    "K": {
        "temperature_range": (3700, 5200),
        "mass_range": (0.45, 0.8),
        "luminosity_range": (0.08, 0.6),
        "color": (1.0, 0.8, 0.6),
        "lifetime_myr": 50000,
        "main_sequence_fraction": 0.121,
        "examples": ["Arcturus", "Aldebaran", "Alpha Centauri B"]
    },
    "M": {
        "temperature_range": (2400, 3700),
        "mass_range": (0.08, 0.45),
        "luminosity_range": (0.0001, 0.08),
        "color": (1.0, 0.6, 0.4),
        "lifetime_myr": 100000,
        "main_sequence_fraction": 0.763,
        "examples": ["Proxima Centauri", "Barnard's Star", "Wolf 359"]
    }
}

# Galaxy Classification (Hubble Sequence)
GALAXY_TYPES = {
    "E0": {
        "description": "Elliptical, spherical",
        "ellipticity": 0.0,
        "bulge_fraction": 1.0,
        "disk_fraction": 0.0,
        "star_formation_rate": 0.1,
        "typical_mass": 1e11,
        "color_index": 0.9,
        "examples": ["M87", "NGC 4472"]
    },
    "E3": {
        "description": "Elliptical, moderate flattening",
        "ellipticity": 0.3,
        "bulge_fraction": 1.0,
        "disk_fraction": 0.0,
        "star_formation_rate": 0.2,
        "typical_mass": 8e10,
        "color_index": 0.85,
        "examples": ["NGC 4374", "NGC 4649"]
    },
    "E7": {
        "description": "Elliptical, highly flattened",
        "ellipticity": 0.7,
        "bulge_fraction": 1.0,
        "disk_fraction": 0.0,
        "star_formation_rate": 0.3,
        "typical_mass": 5e10,
        "color_index": 0.8,
        "examples": ["NGC 3115", "NGC 5866"]
    },
    "S0": {
        "description": "Lenticular, no spiral arms",
        "ellipticity": 0.4,
        "bulge_fraction": 0.6,
        "disk_fraction": 0.4,
        "star_formation_rate": 0.5,
        "typical_mass": 7e10,
        "color_index": 0.75,
        "examples": ["NGC 5866", "Cartwheel Galaxy"]
    },
    "Sa": {
        "description": "Spiral, tightly wound arms",
        "ellipticity": 0.2,
        "bulge_fraction": 0.4,
        "disk_fraction": 0.6,
        "star_formation_rate": 2.0,
        "typical_mass": 1e11,
        "color_index": 0.6,
        "examples": ["M81", "NGC 4594"]
    },
    "Sb": {
        "description": "Spiral, moderately wound arms",
        "ellipticity": 0.3,
        "bulge_fraction": 0.3,
        "disk_fraction": 0.7,
        "star_formation_rate": 3.0,
        "typical_mass": 8e10,
        "color_index": 0.5,
        "examples": ["M31 (Andromeda)", "M101"]
    },
    "Sc": {
        "description": "Spiral, loosely wound arms",
        "ellipticity": 0.4,
        "bulge_fraction": 0.2,
        "disk_fraction": 0.8,
        "star_formation_rate": 5.0,
        "typical_mass": 5e10,
        "color_index": 0.4,
        "examples": ["Milky Way", "M33"]
    },
    "SBa": {
        "description": "Barred spiral, tight arms",
        "ellipticity": 0.2,
        "bulge_fraction": 0.35,
        "disk_fraction": 0.65,
        "star_formation_rate": 2.5,
        "typical_mass": 9e10,
        "color_index": 0.55,
        "examples": ["NGC 1365", "NGC 7479"]
    },
    "SBb": {
        "description": "Barred spiral, moderate arms",
        "ellipticity": 0.3,
        "bulge_fraction": 0.25,
        "disk_fraction": 0.75,
        "star_formation_rate": 4.0,
        "typical_mass": 7e10,
        "color_index": 0.45,
        "examples": ["NGC 1300", "M95"]
    },
    "SBc": {
        "description": "Barred spiral, loose arms",
        "ellipticity": 0.4,
        "bulge_fraction": 0.15,
        "disk_fraction": 0.85,
        "star_formation_rate": 6.0,
        "typical_mass": 4e10,
        "color_index": 0.35,
        "examples": ["NGC 175", "M91"]
    },
    "Irr": {
        "description": "Irregular galaxy",
        "ellipticity": 0.6,
        "bulge_fraction": 0.1,
        "disk_fraction": 0.9,
        "star_formation_rate": 8.0,
        "typical_mass": 1e9,
        "color_index": 0.2,
        "examples": ["Large Magellanic Cloud", "Small Magellanic Cloud", "M82"]
    }
}

# HR Diagram Parameters
HR_DIAGRAM_PARAMS = {
    "main_sequence": {
        "mass_range": (0.08, 120),
        "temperature_range": (2400, 50000),
        "luminosity_range": (1e-4, 1e6),
        "lifetime_range": (1e6, 1e11),  # years
        "color_range": (-0.4, 2.0)  # B-V color index
    },
    "red_giants": {
        "mass_range": (0.5, 8),
        "temperature_range": (3000, 5000),
        "luminosity_range": (10, 1000),
        "lifetime_range": (1e6, 1e9),
        "color_range": (1.0, 2.0)
    },
    "white_dwarfs": {
        "mass_range": (0.17, 1.4),
        "temperature_range": (4000, 150000),
        "luminosity_range": (1e-4, 100),
        "lifetime_range": (1e9, 1e15),
        "color_range": (-0.5, 1.5)
    },
    "supergiants": {
        "mass_range": (8, 120),
        "temperature_range": (3000, 50000),
        "luminosity_range": (1000, 1000000),
        "lifetime_range": (1e6, 5e7),
        "color_range": (0.0, 2.0)
    }
}

# ==========================================
# Data Generation Functions
# ==========================================

def create_sample_stellar_data(n_stars: int = 1000, spectral_classes: Optional[List[str]] = None) -> Dict[str, np.ndarray]:
    """
    Create sample stellar data with realistic distributions.
    
    Args:
        n_stars: Number of stars to generate
        spectral_classes: List of spectral classes to include (default: all)
        
    Returns:
        Dict containing stellar parameters as numpy arrays
    """
    if spectral_classes is None:
        spectral_classes = list(STELLAR_CLASSIFICATION.keys())
    
    # Generate spectral classes based on realistic frequencies
    class_weights = [STELLAR_CLASSIFICATION[sc]["main_sequence_fraction"] for sc in spectral_classes]
    class_weights = np.array(class_weights) / np.sum(class_weights)
    
    stellar_classes = np.random.choice(spectral_classes, size=n_stars, p=class_weights)
    
    # Initialize arrays
    temperatures = np.zeros(n_stars)
    masses = np.zeros(n_stars)
    luminosities = np.zeros(n_stars)
    colors_r = np.zeros(n_stars)
    colors_g = np.zeros(n_stars) 
    colors_b = np.zeros(n_stars)
    lifetimes = np.zeros(n_stars)
    
    # Generate properties for each star
    for i, sc in enumerate(stellar_classes):
        params = STELLAR_CLASSIFICATION[sc]
        
        # Temperature (log-normal distribution within range)
        temp_min, temp_max = params["temperature_range"]
        temperatures[i] = np.random.uniform(temp_min, temp_max)
        
        # Mass (log-normal distribution)
        mass_min, mass_max = params["mass_range"]
        masses[i] = np.random.lognormal(
            np.log(np.sqrt(mass_min * mass_max)), 
            0.3
        )
        masses[i] = np.clip(masses[i], mass_min, mass_max)
        
        # Luminosity (mass-luminosity relation with scatter)
        if masses[i] < 0.43:
            lum_base = 0.23 * (masses[i] ** 2.3)
        elif masses[i] < 2:
            lum_base = masses[i] ** 4
        elif masses[i] < 20:
            lum_base = 1.5 * (masses[i] ** 3.5)
        else:
            lum_base = 32000 * masses[i]
            
        luminosities[i] = lum_base * np.random.lognormal(0, 0.2)
        
        # Color
        colors_r[i], colors_g[i], colors_b[i] = params["color"]
        
        # Add scatter to colors
        color_scatter = 0.1
        colors_r[i] += np.random.normal(0, color_scatter)
        colors_g[i] += np.random.normal(0, color_scatter)
        colors_b[i] += np.random.normal(0, color_scatter)
        
        # Lifetime
        lifetimes[i] = params["lifetime_myr"] * np.random.lognormal(0, 0.3)
    
    # Return as numpy arrays
    data = {
        "spectral_class": stellar_classes,
        "temperature": temperatures.astype(np.float32),
        "mass": masses.astype(np.float32),
        "luminosity": luminosities.astype(np.float32),
        "color_r": colors_r.astype(np.float32),
        "color_g": colors_g.astype(np.float32),
        "color_b": colors_b.astype(np.float32),
        "lifetime_myr": lifetimes.astype(np.float32),
    }
    
    return data

def create_sample_hr_diagram_data(n_stars: int = 2000) -> Dict[str, np.ndarray]:
    """
    Create sample HR diagram data with all stellar evolution phases.
    
    Args:
        n_stars: Number of stars to generate
        
    Returns:
        Dict containing HR diagram data as numpy arrays
    """
    # Distribution of stellar evolution phases
    phase_weights = {
        "main_sequence": 0.85,
        "red_giants": 0.08,
        "white_dwarfs": 0.05,
        "supergiants": 0.02
    }
    
    phases = np.random.choice(
        list(phase_weights.keys()), 
        size=n_stars, 
        p=list(phase_weights.values())
    )
    
    # Initialize arrays
    temperatures = np.zeros(n_stars)
    luminosities = np.zeros(n_stars)
    masses = np.zeros(n_stars)
    colors = np.zeros(n_stars)
    
    # Generate properties for each stellar phase
    for i, phase in enumerate(phases):
        params = HR_DIAGRAM_PARAMS[phase]
        
        # Temperature (log-uniform in range)
        temp_min, temp_max = params["temperature_range"]
        temperatures[i] = np.random.uniform(temp_min, temp_max)
        
        # Luminosity (log-uniform in range)
        lum_min, lum_max = params["luminosity_range"]
        luminosities[i] = np.random.uniform(
            np.log10(lum_min), np.log10(lum_max)
        )
        luminosities[i] = 10 ** luminosities[i]
        
        # Mass (log-uniform in range)
        mass_min, mass_max = params["mass_range"]
        masses[i] = np.random.uniform(mass_min, mass_max)
        
        # Color index (B-V)
        color_min, color_max = params["color_range"]
        colors[i] = np.random.uniform(color_min, color_max)
    
    # Return as numpy arrays
    data = {
        "evolution_phase": phases,
        "temperature": temperatures.astype(np.float32),
        "luminosity": luminosities.astype(np.float32),
        "mass": masses.astype(np.float32),
        "color_bv": colors.astype(np.float32),
        "log_temperature": np.log10(temperatures).astype(np.float32),
        "log_luminosity": np.log10(luminosities).astype(np.float32),
    }
    
    return data

# ==========================================
# Utility Functions
# ==========================================

def get_stellar_data(spectral_class: str) -> Dict[str, Any]:
    """Get stellar data for a specific spectral class."""
    return STELLAR_CLASSIFICATION.get(spectral_class, STELLAR_CLASSIFICATION["G"])

def get_galaxy_config(galaxy_type: str) -> Dict[str, Any]:
    """Get galaxy configuration for a specific type."""
    return GALAXY_TYPES.get(galaxy_type, GALAXY_TYPES["Sb"])

def validate_hr_diagram_data(data: Dict[str, np.ndarray]) -> bool:
    """
    Validate HR diagram data for physical consistency.
    
    Args:
        data: HR diagram data dict with numpy arrays
        
    Returns:
        True if data is valid
    """
    required_keys = ["temperature", "luminosity"]
    
    # Check required keys
    for key in required_keys:
        if key not in data:
            print(f"Missing required key: {key}")
            return False
    
    # Check data ranges
    temp = data["temperature"]
    lum = data["luminosity"]
    
    if np.any(temp < 1000) or np.any(temp > 100000):
        print("Warning: Temperature values outside expected range (1000-100000 K)")
    
    if np.any(lum < 1e-6) or np.any(lum > 1e7):
        print("Warning: Luminosity values outside expected range (1e-6 - 1e7 solar)")
    
    # Check for NaN or infinite values
    if np.any(np.isnan(temp)) or np.any(np.isinf(temp)):
        print("Error: Invalid temperature values (NaN or Inf)")
        return False
        
    if np.any(np.isnan(lum)) or np.any(np.isinf(lum)):
        print("Error: Invalid luminosity values (NaN or Inf)")
        return False
    
    print("HR diagram data validation passed")
    return True

def get_stellar_color_from_temperature(temperature: Union[float, np.ndarray]) -> Union[tuple, np.ndarray]:
    """
    Convert stellar temperature to RGB color using blackbody approximation
    
    Args:
        temperature: Temperature in Kelvin (single value or array)
    
    Returns:
        RGB color(s) as tuple or array
    """
    # Simplified blackbody color conversion
    is_array = isinstance(temperature, np.ndarray)
    
    if not is_array:
        temperature = np.array([temperature])
    
    # Initialize RGB arrays
    red = np.ones_like(temperature, dtype=np.float32)
    green = np.ones_like(temperature, dtype=np.float32)
    blue = np.ones_like(temperature, dtype=np.float32)
    
    # Red calculation
    mask = temperature < 6600
    red[mask] = 255
    red[~mask] = 329.698727446 * ((temperature[~mask] / 100 - 60) ** -0.1332047592)
    
    # Green calculation
    mask = temperature < 6600
    green[mask] = 99.4708025861 * np.log(temperature[mask] / 100) - 161.1195681661
    green[~mask] = 288.1221695283 * ((temperature[~mask] / 100 - 60) ** -0.0755148492)
    
    # Blue calculation
    mask = temperature >= 6600
    blue[mask] = 255
    mask = temperature < 1900
    blue[mask] = 0
    mask = (temperature >= 1900) & (temperature < 6600)
    blue[mask] = 138.5177312231 * np.log(temperature[mask] / 100 - 10) - 305.0447927307
    
    # Normalize to 0-1 range
    red = np.clip(red / 255.0, 0, 1)
    green = np.clip(green / 255.0, 0, 1)
    blue = np.clip(blue / 255.0, 0, 1)
    
    if is_array:
        return np.column_stack([red, green, blue])
    else:
        return (float(red[0]), float(green[0]), float(blue[0]))

def calculate_absolute_magnitude(luminosity: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Calculate absolute magnitude from luminosity (in solar units)
    
    Args:
        luminosity: Luminosity in solar units
    
    Returns:
        Absolute magnitude
    """
    # Absolute magnitude of the Sun
    solar_abs_mag = 4.83
    
    # Magnitude formula: M = M_sun - 2.5 * log10(L/L_sun)
    return solar_abs_mag - 2.5 * np.log10(luminosity)

def calculate_distance_modulus(apparent_mag: Union[float, np.ndarray], 
                              absolute_mag: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Calculate distance modulus from apparent and absolute magnitudes
    
    Args:
        apparent_mag: Apparent magnitude
        absolute_mag: Absolute magnitude
    
    Returns:
        Distance modulus (apparent - absolute)
    """
    return apparent_mag - absolute_mag

def distance_modulus_to_parsecs(distance_modulus: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert distance modulus to distance in parsecs
    
    Args:
        distance_modulus: Distance modulus
    
    Returns:
        Distance in parsecs
    """
    return 10 ** ((distance_modulus + 5) / 5)
