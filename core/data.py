"""
Data Access for Astronomical Extension
=====================================

Data loading and processing for astronomical visualization.
No external dependencies - self-contained implementation.
"""

import numpy as np
from typing import Dict, Any, Tuple, List, Optional

# ==========================================
# Sample Data Generation
# ==========================================

def list_available_surveys():
    """
    List all available survey names for demonstration
    
    Returns:
        list: Available survey identifiers
    """
    return ["gaia_demo", "sdss_sample", "nsa_preview", "synthetic_stars", "synthetic_galaxies"]

def load_survey_data(survey: str, max_samples: int = 10000) -> Dict[str, Any]:
    """
    Load sample astronomical data for visualization
    
    Args:
        survey: Survey identifier
        max_samples: Maximum number of objects to generate
    
    Returns:
        dict: Sample data with coordinates and properties
    """
    if survey == "gaia_demo":
        return _generate_gaia_sample(max_samples)
    elif survey == "sdss_sample":
        return _generate_sdss_sample(max_samples)
    elif survey == "synthetic_stars":
        return _generate_stellar_sample(max_samples)
    elif survey == "synthetic_galaxies":
        return _generate_galaxy_sample(max_samples)
    else:
        # Default to stellar sample
        return _generate_stellar_sample(max_samples)

def get_coordinates_and_features(survey: str, max_samples: int = 10000) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Load survey data and extract coordinates and features for visualization
    
    Args:
        survey: Survey identifier
        max_samples: Maximum number of objects
    
    Returns:
        tuple: (coordinates_array, features_dict)
            coordinates: numpy array of shape (N, 3) - x, y, z positions
            features: dict with magnitudes, colors, types, etc.
    """
    data = load_survey_data(survey, max_samples)
    
    # Extract coordinates
    coords = np.column_stack([
        data['x'],
        data['y'], 
        data['z']
    ])
    
    # Extract features
    features = {key: value for key, value in data.items() if key not in ['x', 'y', 'z']}
    
    return coords, features

# ==========================================
# Sample Data Generators
# ==========================================

def _generate_gaia_sample(count: int) -> Dict[str, Any]:
    """Generate Gaia-like stellar data"""
    np.random.seed(42)  # Reproducible results
    
    # Galactic disk distribution
    radius = np.random.exponential(3.0, count) * 1000  # pc
    theta = np.random.uniform(0, 2*np.pi, count)
    z_height = np.random.normal(0, 200, count)  # pc
    
    # Cartesian coordinates
    x = radius * np.cos(theta)
    y = radius * np.sin(theta) 
    z = z_height
    
    # Stellar properties
    magnitude = np.random.uniform(5, 20, count)  # G magnitude
    bp_rp = np.random.normal(0.8, 0.5, count)  # Color index
    parallax = 1000 / (radius + 1)  # mas (distance modulus)
    
    return {
        'x': x,
        'y': y,
        'z': z,
        'magnitude': magnitude,
        'bp_rp_color': bp_rp,
        'parallax': parallax,
        'survey_type': 'gaia'
    }

def _generate_sdss_sample(count: int) -> Dict[str, Any]:
    """Generate SDSS-like galaxy data"""
    np.random.seed(123)
    
    # Large scale structure distribution
    x = np.random.normal(0, 50000, count)  # Mpc
    y = np.random.normal(0, 50000, count) 
    z = np.random.uniform(0, 200000, count)  # Redshift space
    
    # Galaxy properties
    r_magnitude = np.random.uniform(14, 22, count)
    g_r_color = np.random.normal(0.6, 0.3, count)
    redshift = z / 300000  # Simple redshift approximation
    
    # Galaxy types
    galaxy_types = np.random.choice(['spiral', 'elliptical', 'irregular'], count, 
                                   p=[0.6, 0.3, 0.1])
    
    return {
        'x': x,
        'y': y,
        'z': z,
        'r_magnitude': r_magnitude,
        'g_r_color': g_r_color,
        'redshift': redshift,
        'galaxy_type': galaxy_types,
        'survey_type': 'sdss'
    }

def _generate_stellar_sample(count: int) -> Dict[str, Any]:
    """Generate synthetic stellar population"""
    np.random.seed(456)
    
    # Spherical distribution
    radius = np.random.uniform(1, 100, count)
    theta = np.random.uniform(0, 2*np.pi, count)
    phi = np.random.uniform(0, np.pi, count)
    
    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)
    
    # Stellar classification
    spectral_types = np.random.choice(['O', 'B', 'A', 'F', 'G', 'K', 'M'], count,
                                     p=[0.01, 0.05, 0.15, 0.2, 0.25, 0.25, 0.09])
    
    # Magnitude based on spectral type
    magnitude_map = {'O': 3, 'B': 4, 'A': 5, 'F': 6, 'G': 7, 'K': 8, 'M': 10}
    magnitude = np.array([magnitude_map[st] + np.random.normal(0, 1) for st in spectral_types])
    
    # Color based on spectral type  
    color_map = {'O': -0.3, 'B': 0.0, 'A': 0.2, 'F': 0.4, 'G': 0.6, 'K': 1.0, 'M': 1.5}
    color = np.array([color_map[st] + np.random.normal(0, 0.1) for st in spectral_types])
    
    return {
        'x': x,
        'y': y,
        'z': z,
        'magnitude': magnitude,
        'color': color,
        'spectral_type': spectral_types,
        'survey_type': 'synthetic_stars'
    }

def _generate_galaxy_sample(count: int) -> Dict[str, Any]:
    """Generate synthetic galaxy distribution"""
    np.random.seed(789)
    
    # Cosmic web-like distribution
    # Create filamentary structure
    filament_count = max(1, count // 100)
    galaxies_per_filament = count // filament_count
    
    x_coords = []
    y_coords = []
    z_coords = []
    
    for i in range(filament_count):
        # Random filament direction
        direction = np.random.uniform(-1, 1, 3)
        direction = direction / np.linalg.norm(direction)
        
        # Filament length
        length = np.random.uniform(10000, 50000)  # Mpc
        
        # Generate galaxies along filament
        t_values = np.random.uniform(0, length, galaxies_per_filament)
        
        # Add scatter perpendicular to filament
        scatter = np.random.normal(0, 1000, (galaxies_per_filament, 3))  # Mpc
        
        # Calculate positions
        for j, t in enumerate(t_values):
            position = t * direction + scatter[j]
            x_coords.append(position[0])
            y_coords.append(position[1])
            z_coords.append(position[2])
    
    # Trim to exact count
    x_coords = np.array(x_coords[:count])
    y_coords = np.array(y_coords[:count])
    z_coords = np.array(z_coords[:count])
    
    # Galaxy properties
    magnitude = np.random.uniform(12, 20, count)
    color = np.random.normal(0.5, 0.3, count)
    mass = np.random.lognormal(10, 1, count)  # Solar masses
    
    return {
        'x': x_coords,
        'y': y_coords,
        'z': z_coords,
        'magnitude': magnitude,
        'color': color,
        'stellar_mass': mass,
        'survey_type': 'synthetic_galaxies'
    }

# ==========================================
# Data Utility Functions
# ==========================================

def get_data_statistics(survey: str, max_samples: int = 10000) -> Dict[str, Any]:
    """Get statistical summary of survey data"""
    data = load_survey_data(survey, max_samples)
    coords, features = get_coordinates_and_features(survey, max_samples)
    
    stats = {
        'count': len(coords),
        'spatial_extent': {
            'x_range': (float(np.min(coords[:, 0])), float(np.max(coords[:, 0]))),
            'y_range': (float(np.min(coords[:, 1])), float(np.max(coords[:, 1]))),
            'z_range': (float(np.min(coords[:, 2])), float(np.max(coords[:, 2]))),
        },
        'survey_type': data.get('survey_type', 'unknown'),
        'available_features': list(features.keys())
    }
    
    return stats

def validate_survey_data(survey: str) -> Dict[str, Any]:
    """Validate that survey data can be loaded properly"""
    try:
        coords, features = get_coordinates_and_features(survey, 100)  # Small test sample
        
        return {
            'valid': True,
            'coordinates_shape': coords.shape,
            'feature_count': len(features),
            'coordinate_range': {
                'x': (float(np.min(coords[:, 0])), float(np.max(coords[:, 0]))),
                'y': (float(np.min(coords[:, 1])), float(np.max(coords[:, 1]))),
                'z': (float(np.min(coords[:, 2])), float(np.max(coords[:, 2]))),
            }
        }
    except Exception as e:
        return {
            'valid': False,
            'error': str(e)
        }

def create_sample_dataset(survey_type: str = "mixed", count: int = 1000) -> Dict[str, Any]:
    """Create a mixed sample dataset for testing"""
    if survey_type == "mixed":
        # Combine different object types
        star_data = _generate_stellar_sample(count // 2)
        galaxy_data = _generate_galaxy_sample(count // 2)
        
        # Combine data
        combined_data = {
            'x': np.concatenate([star_data['x'], galaxy_data['x']]),
            'y': np.concatenate([star_data['y'], galaxy_data['y']]),
            'z': np.concatenate([star_data['z'], galaxy_data['z']]),
            'magnitude': np.concatenate([star_data['magnitude'], galaxy_data['magnitude']]),
            'object_type': (['star'] * len(star_data['x'])) + (['galaxy'] * len(galaxy_data['x']))
        }
        
        return combined_data
    else:
        return load_survey_data(survey_type, count)
