"""
FITS Utilities for AstroPhot Integration
=======================================

Enhanced FITS utilities with complete Astropy integration.
Supports WCS, coordinate transformations, I/O operations, and astronomical calculations.
"""

import numpy as np
from typing import Dict, List, Any, Optional, Union, Tuple
from pathlib import Path

# Astropy imports - complete integration following standards
import astropy.units as u
import astropy.constants as const
import astropy.coordinates as coord
import astropy.time as time
import astropy.io.fits as fits
import astropy.wcs as wcs
import astropy.table as table
from astropy.coordinates import SkyCoord
from astropy.time import Time

# Import our enhanced Astropy integration
from .astropy_integration import (
    AstronomicalConstants, CoordinateSystem, AstronomicalTime,
    WorldCoordinateSystem, AstropyIO, UnitConversions, TableOperations
)

# FITS file extensions
FITS_EXTENSIONS = ['.fits', '.fit', '.fts']

def validate_fits_file(filepath: Union[str, Path]) -> bool:
    """
    Enhanced FITS file validation with Astropy verification.
    
    Args:
        filepath: Path to the file
    
    Returns:
        True if valid FITS file
    """
    filepath = Path(filepath)
    
    # Check extension
    if filepath.suffix.lower() not in FITS_EXTENSIONS:
        return False
    
    # Check if file exists
    if not filepath.exists():
        return False
    
    # Check file size (FITS files should be at least a few KB)
    if filepath.stat().st_size < 1024:
        return False
    
    # Validate with Astropy
    try:
        with fits.open(filepath) as hdul:
            # Check if it's a valid FITS file
            if len(hdul) == 0:
                return False
            
            # Check primary HDU
            primary = hdul[0]
            if primary.header is None:
                return False
            
            return True
    except Exception:
        return False

def get_fits_info(data: np.ndarray, header: Optional[fits.Header] = None, 
                 filepath: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """
    Extract comprehensive information from FITS data with Astropy integration.
    
    Args:
        data: FITS data as numpy array
        header: FITS header (optional)
        filepath: Path to FITS file (optional)
    
    Returns:
        Dictionary with comprehensive FITS information
    """
    info = {
        'shape': data.shape,
        'dtype': str(data.dtype),
        'size_mb': data.nbytes / (1024 * 1024),
        'dimensions': len(data.shape),
    }
    
    # Basic statistics (handle NaN values)
    valid_data = data[~np.isnan(data)] if np.any(np.isnan(data)) else data
    
    if len(valid_data) > 0:
        info.update({
            'min_value': float(np.min(valid_data)),
            'max_value': float(np.max(valid_data)),
            'mean_value': float(np.mean(valid_data)),
            'std_value': float(np.std(valid_data)),
            'nan_pixels': int(np.sum(np.isnan(data))),
            'total_pixels': int(data.size),
        })
    
    # Enhanced header analysis
    if header is not None:
        info.update(analyze_fits_header(header))
    
    # Load header from file if provided
    if filepath is not None and header is None:
        try:
            header = AstropyIO.read_fits_header(filepath)
            info.update(analyze_fits_header(header))
        except Exception as e:
            info['header_error'] = str(e)
    
    # WCS information
    if header is not None:
        try:
            wcs_info = analyze_wcs(header)
            info['wcs'] = wcs_info
        except Exception as e:
            info['wcs_error'] = str(e)
    
    # Suggest visualization parameters
    info['visualization'] = suggest_visualization_params(data)
    
    # Object type estimation
    info['object_type'] = estimate_object_type(data)
    
    return info

def analyze_fits_header(header: fits.Header) -> Dict[str, Any]:
    """
    Comprehensive FITS header analysis with astronomical interpretation.
    
    Args:
        header: FITS header
    
    Returns:
        Dictionary with header analysis
    """
    analysis = {}
    
    # Basic header info
    analysis['n_keywords'] = len(header)
    analysis['instrument'] = header.get('INSTRUME', 'Unknown')
    analysis['telescope'] = header.get('TELESCOP', 'Unknown')
    analysis['observer'] = header.get('OBSERVER', 'Unknown')
    analysis['object_name'] = header.get('OBJECT', 'Unknown')
    
    # Observational information
    analysis['date_obs'] = header.get('DATE-OBS', 'Unknown')
    analysis['exposure_time'] = header.get('EXPTIME', header.get('EXPOSURE', None))
    analysis['filter'] = header.get('FILTER', header.get('FILTNAM', 'Unknown'))
    
    # Coordinate information
    analysis['ra'] = header.get('RA', header.get('CRVAL1', None))
    analysis['dec'] = header.get('DEC', header.get('CRVAL2', None))
    analysis['equinox'] = header.get('EQUINOX', header.get('EPOCH', None))
    
    # Image characteristics
    analysis['naxis'] = header.get('NAXIS', 0)
    if analysis['naxis'] >= 2:
        analysis['naxis1'] = header.get('NAXIS1', 0)
        analysis['naxis2'] = header.get('NAXIS2', 0)
    
    # Photometric information
    analysis['gain'] = header.get('GAIN', None)
    analysis['readnoise'] = header.get('RDNOISE', header.get('READNOIS', None))
    analysis['airmass'] = header.get('AIRMASS', None)
    analysis['seeing'] = header.get('SEEING', header.get('FWHM', None))
    
    # Data scaling
    analysis['bzero'] = header.get('BZERO', 0.0)
    analysis['bscale'] = header.get('BSCALE', 1.0)
    analysis['bunit'] = header.get('BUNIT', 'Unknown')
    
    # Time information with Astropy
    if analysis['date_obs'] != 'Unknown':
        try:
            obs_time = AstronomicalTime.create_time(analysis['date_obs'])
            analysis['mjd'] = obs_time.mjd
            analysis['jd'] = obs_time.jd
            analysis['time_scale'] = obs_time.scale
        except Exception:
            pass
    
    # Coordinate object if available
    if analysis['ra'] is not None and analysis['dec'] is not None:
        try:
            coord_obj = CoordinateSystem.create_skycoord(
                analysis['ra'], analysis['dec'], unit='deg'
            )
            analysis['coordinates'] = {
                'ra_deg': coord_obj.ra.deg,
                'dec_deg': coord_obj.dec.deg,
                'ra_hms': coord_obj.ra.to_string(unit=u.hour),
                'dec_dms': coord_obj.dec.to_string(unit=u.deg),
                'galactic_l': coord_obj.galactic.l.deg,
                'galactic_b': coord_obj.galactic.b.deg,
            }
        except Exception:
            pass
    
    return analysis

def analyze_wcs(header: fits.Header) -> Dict[str, Any]:
    """
    Comprehensive WCS analysis using Astropy.
    
    Args:
        header: FITS header with WCS information
    
    Returns:
        Dictionary with WCS analysis
    """
    try:
        wcs_obj = WorldCoordinateSystem.create_wcs_from_fits(header)
        
        analysis = {
            'has_wcs': True,
            'wcs_type': 'celestial' if wcs_obj.has_celestial else 'other',
            'coordinate_system': wcs_obj.wcs.ctype,
            'reference_pixel': wcs_obj.wcs.crpix,
            'reference_coordinates': wcs_obj.wcs.crval,
            'pixel_scale': wcs_obj.wcs.cdelt,
            'rotation_matrix': wcs_obj.wcs.cd if hasattr(wcs_obj.wcs, 'cd') else None,
        }
        
        # Calculate pixel scales
        if wcs_obj.has_celestial:
            try:
                pixel_scales = WorldCoordinateSystem.calculate_pixel_scale(wcs_obj)
                analysis['pixel_scale_arcsec'] = pixel_scales
                analysis['field_of_view_arcmin'] = [
                    header.get('NAXIS1', 0) * pixel_scales[0] / 60,
                    header.get('NAXIS2', 0) * pixel_scales[1] / 60
                ]
            except Exception:
                pass
        
        # Projection information
        analysis['projection'] = wcs_obj.wcs.ctype[0].split('-')[-1] if len(wcs_obj.wcs.ctype) > 0 else 'Unknown'
        
        return analysis
        
    except Exception as e:
        return {
            'has_wcs': False,
            'error': str(e)
        }

def suggest_visualization_params(data: np.ndarray) -> Dict[str, Any]:
    """
    Enhanced visualization parameter suggestions with astronomical context.
    
    Args:
        data: FITS data array
    
    Returns:
        Suggested visualization parameters
    """
    valid_data = data[~np.isnan(data)] if np.any(np.isnan(data)) else data
    
    if len(valid_data) == 0:
        return {'error': 'No valid data'}
    
    # Calculate percentiles for better visualization
    percentiles = [0.1, 1, 5, 25, 50, 75, 95, 99, 99.9]
    p_values = np.percentile(valid_data, percentiles)
    
    # Determine if data is mostly positive (typical for astronomical images)
    positive_fraction = np.sum(valid_data > 0) / len(valid_data)
    negative_fraction = np.sum(valid_data < 0) / len(valid_data)
    
    # Signal-to-noise estimation
    background_level = p_values[3]  # 25th percentile as background estimate
    noise_level = np.std(valid_data[valid_data < p_values[6]])  # Std of lower values
    snr = (p_values[7] - background_level) / noise_level if noise_level > 0 else 1.0
    
    suggestions = {
        'data_range': (float(np.min(valid_data)), float(np.max(valid_data))),
        'suggested_range': (float(p_values[1]), float(p_values[7])),  # 1% to 99%
        'conservative_range': (float(p_values[2]), float(p_values[6])),  # 5% to 95%
        'robust_range': (float(p_values[0]), float(p_values[8])),  # 0.1% to 99.9%
        'positive_fraction': float(positive_fraction),
        'negative_fraction': float(negative_fraction),
        'snr_estimate': float(snr),
        'background_level': float(background_level),
        'noise_level': float(noise_level),
        'percentiles': dict(zip(percentiles, p_values.astype(float))),
    }
    
    # Suggest color mapping based on data characteristics
    if positive_fraction > 0.95:
        suggestions['suggested_colormap'] = 'BLACKBODY'
        suggestions['reason'] = 'Mostly positive values - likely emission/intensity data'
    elif negative_fraction > 0.3:
        suggestions['suggested_colormap'] = 'RAINBOW'
        suggestions['reason'] = 'Significant negative values - likely difference/residual image'
    elif snr > 10:
        suggestions['suggested_colormap'] = 'VIRIDIS'
        suggestions['reason'] = 'High SNR data - scientific colormap recommended'
    else:
        suggestions['suggested_colormap'] = 'GRAYSCALE'
        suggestions['reason'] = 'Low SNR data - grayscale for better contrast'
    
    # Suggest displacement scale based on data characteristics
    data_range = p_values[7] - p_values[1]  # 99% - 1% range
    if data_range > 0:
        # Scale displacement to reasonable range (0.1 to 3.0 Blender units)
        suggestions['displacement_scale'] = float(np.clip(2.0 / (data_range / 1000), 0.1, 3.0))
    else:
        suggestions['displacement_scale'] = 1.0
    
    return suggestions

def estimate_object_type(data: np.ndarray) -> Dict[str, Any]:
    """
    Enhanced astronomical object type estimation.
    
    Args:
        data: FITS data array
    
    Returns:
        Dictionary with estimated object type and confidence
    """
    valid_data = data[~np.isnan(data)] if np.any(np.isnan(data)) else data
    
    if len(valid_data) == 0:
        return {'type': 'unknown', 'confidence': 0.0, 'reason': 'No valid data'}
    
    # Calculate enhanced statistics
    mean_val = np.mean(valid_data)
    std_val = np.std(valid_data)
    max_val = np.max(valid_data)
    
    # Calculate signal-to-noise ratio
    snr = mean_val / std_val if std_val > 0 else 0
    
    # Analyze data distribution
    positive_fraction = np.sum(valid_data > mean_val) / len(valid_data)
    
    # Enhanced classification logic
    if snr > 50 and max_val / mean_val > 100:
        return {
            'type': 'stellar_field_bright',
            'confidence': 0.9,
            'reason': 'Very high SNR with bright point sources',
            'snr': float(snr),
            'suggested_analysis': ['PSF photometry', 'Stellar classification']
        }
    elif snr > 20:
        return {
            'type': 'stellar_field',
            'confidence': 0.8,
            'reason': 'High SNR with point sources',
            'snr': float(snr),
            'suggested_analysis': ['Aperture photometry', 'Source detection']
        }
    elif snr > 10 and positive_fraction > 0.7:
        return {
            'type': 'galaxy',
            'confidence': 0.7,
            'reason': 'Moderate SNR with extended structure',
            'snr': float(snr),
            'suggested_analysis': ['Surface brightness profile', 'Morphological analysis']
        }
    elif snr > 5:
        return {
            'type': 'nebula',
            'confidence': 0.6,
            'reason': 'Extended emission structure',
            'snr': float(snr),
            'suggested_analysis': ['Emission line analysis', 'Contour mapping']
        }
    elif abs(mean_val) < std_val:
        return {
            'type': 'noise_or_background',
            'confidence': 0.9,
            'reason': 'Low signal relative to noise',
            'snr': float(snr),
            'suggested_analysis': ['Bias/dark correction', 'Flat fielding']
        }
    else:
        return {
            'type': 'unknown_extended',
            'confidence': 0.4,
            'reason': 'Unclear data characteristics',
            'snr': float(snr),
            'suggested_analysis': ['Visual inspection', 'Multi-wavelength comparison']
        }

def prepare_fits_for_blender(data: np.ndarray, normalize: bool = True, 
                           clip_percentiles: Optional[Tuple[float, float]] = None,
                           scaling_method: str = 'linear') -> np.ndarray:
    """
    Enhanced FITS data preparation for Blender visualization.
    
    Args:
        data: FITS data array
        normalize: Whether to normalize to 0-1 range
        clip_percentiles: Tuple of (low, high) percentiles for clipping
        scaling_method: Scaling method ('linear', 'sqrt', 'log', 'asinh')
    
    Returns:
        Processed data ready for Blender
    """
    processed = data.copy()
    
    # Handle NaN values by replacing with minimum
    if np.any(np.isnan(processed)):
        valid_min = np.nanmin(processed)
        processed[np.isnan(processed)] = valid_min
    
    # Clip to percentiles if specified
    if clip_percentiles is not None:
        low_p, high_p = clip_percentiles
        low_val, high_val = np.percentile(processed, [low_p, high_p])
        processed = np.clip(processed, low_val, high_val)
    
    # Apply scaling method
    if scaling_method == 'sqrt':
        # Square root scaling (good for nebulae)
        processed = processed - np.min(processed)  # Ensure positive
        processed = np.sqrt(processed)
    elif scaling_method == 'log':
        # Logarithmic scaling (good for high dynamic range)
        processed = processed - np.min(processed) + 1  # Ensure positive
        processed = np.log10(processed)
    elif scaling_method == 'asinh':
        # Arcsinh scaling (Lupton et al. - good for mixed signals)
        b = np.percentile(processed, 5)  # Softening parameter
        processed = np.arcsinh(processed / b)
    # 'linear' scaling - no transformation needed
    
    # Normalize to 0-1 range
    if normalize:
        data_min = np.min(processed)
        data_max = np.max(processed)
        
        if data_max > data_min:
            processed = (processed - data_min) / (data_max - data_min)
        else:
            processed = np.zeros_like(processed)
    
    return processed.astype(np.float32)

def generate_astronomical_coordinates(shape: Tuple[int, int], 
                                    center_ra_deg: float = 0.0,
                                    center_dec_deg: float = 0.0,
                                    pixel_scale_arcsec: float = 1.0) -> Dict[str, np.ndarray]:
    """
    Generate astronomical coordinate arrays for FITS data
    
    Args:
        shape: (height, width) of the image
        center_ra_deg: Right ascension of center in degrees
        center_dec_deg: Declination of center in degrees
        pixel_scale_arcsec: Pixel scale in arcseconds per pixel
    
    Returns:
        Dictionary with RA and Dec coordinate arrays
    """
    height, width = shape
    
    # Create pixel coordinate arrays
    x_pixels = np.arange(width) - width / 2
    y_pixels = np.arange(height) - height / 2
    
    # Convert to world coordinates
    pixel_scale_deg = pixel_scale_arcsec / 3600.0  # Convert to degrees
    
    # Create coordinate grids
    x_grid, y_grid = np.meshgrid(x_pixels, y_pixels)
    
    # Convert to RA/Dec (simplified projection)
    ra_array = center_ra_deg + x_grid * pixel_scale_deg
    dec_array = center_dec_deg + y_grid * pixel_scale_deg
    
    return {
        'ra_degrees': ra_array.astype(np.float32),
        'dec_degrees': dec_array.astype(np.float32),
        'x_pixels': x_grid.astype(np.float32),
        'y_pixels': y_grid.astype(np.float32),
        'pixel_scale_arcsec': pixel_scale_arcsec,
        'center_ra_deg': center_ra_deg,
        'center_dec_deg': center_dec_deg,
    }

def create_fits_metadata(filepath: Union[str, Path], data: np.ndarray) -> Dict[str, Any]:
    """
    Create comprehensive metadata for imported FITS file
    
    Args:
        filepath: Path to the FITS file
        data: FITS data array
    
    Returns:
        Complete metadata dictionary
    """
    filepath = Path(filepath)
    
    metadata = {
        'file_info': {
            'filename': filepath.name,
            'filepath': str(filepath),
            'file_size_mb': filepath.stat().st_size / (1024 * 1024),
        },
        'data_info': get_fits_info(data),
        'object_type': estimate_object_type(data),
        'import_settings': {
            'timestamp': str(np.datetime64('now')),
            'blender_units': True,
            'coordinate_system': 'pixel',  # Could be extended to WCS
        }
    }
    
    # Add processing recommendations
    suggestions = suggest_visualization_params(data)
    metadata['visualization_suggestions'] = suggestions
    
    return metadata

def calculate_surface_brightness(data: np.ndarray, pixel_scale_arcsec: float = 1.0) -> np.ndarray:
    """
    Calculate surface brightness in astronomical units
    
    Args:
        data: FITS data array (assumed to be in appropriate flux units)
        pixel_scale_arcsec: Pixel scale in arcseconds per pixel
    
    Returns:
        Surface brightness array
    """
    # Convert to surface brightness per square arcsecond
    pixel_area_arcsec2 = pixel_scale_arcsec ** 2
    surface_brightness = data / pixel_area_arcsec2
    
    # Convert to magnitudes per square arcsecond (if positive flux)
    with np.errstate(divide='ignore', invalid='ignore'):
        mag_per_arcsec2 = -2.5 * np.log10(surface_brightness)
    
    # Handle infinite/invalid values
    mag_per_arcsec2[~np.isfinite(mag_per_arcsec2)] = 30.0  # Faint default
    
    return mag_per_arcsec2.astype(np.float32)

def create_false_color_rgb(data_r: np.ndarray, data_g: np.ndarray, data_b: np.ndarray,
                          scales: Tuple[float, float, float] = (1.0, 1.0, 1.0)) -> np.ndarray:
    """
    Create false-color RGB image from multiple FITS arrays
    
    Args:
        data_r: Red channel data
        data_g: Green channel data  
        data_b: Blue channel data
        scales: Scaling factors for each channel
    
    Returns:
        RGB array ready for Blender
    """
    # Ensure all arrays have the same shape
    shape = data_r.shape
    assert data_g.shape == shape and data_b.shape == shape, "All channels must have same shape"
    
    # Apply scaling
    r_scaled = data_r * scales[0]
    g_scaled = data_g * scales[1]
    b_scaled = data_b * scales[2]
    
    # Normalize each channel independently
    def normalize_channel(channel):
        valid_data = channel[~np.isnan(channel)]
        if len(valid_data) > 0:
            p1, p99 = np.percentile(valid_data, [1, 99])
            normalized = np.clip((channel - p1) / (p99 - p1), 0, 1)
        else:
            normalized = np.zeros_like(channel)
        return normalized
    
    r_norm = normalize_channel(r_scaled)
    g_norm = normalize_channel(g_scaled)
    b_norm = normalize_channel(b_scaled)
    
    # Stack into RGB array
    rgb = np.stack([r_norm, g_norm, b_norm], axis=-1)
    
    return rgb.astype(np.float32)

# Constants for common astronomical filters
COMMON_FILTERS = {
    'U': {'central_wavelength_nm': 365, 'color_rgb': (0.4, 0.0, 1.0)},
    'B': {'central_wavelength_nm': 445, 'color_rgb': (0.0, 0.4, 1.0)},
    'V': {'central_wavelength_nm': 551, 'color_rgb': (0.0, 1.0, 0.0)},
    'R': {'central_wavelength_nm': 658, 'color_rgb': (1.0, 0.4, 0.0)},
    'I': {'central_wavelength_nm': 806, 'color_rgb': (1.0, 0.0, 0.0)},
    'g': {'central_wavelength_nm': 477, 'color_rgb': (0.0, 0.6, 1.0)},  # SDSS
    'r': {'central_wavelength_nm': 623, 'color_rgb': (1.0, 0.6, 0.0)},  # SDSS
    'i': {'central_wavelength_nm': 763, 'color_rgb': (1.0, 0.2, 0.0)},  # SDSS
    'Ha': {'central_wavelength_nm': 656, 'color_rgb': (1.0, 0.0, 0.0)}, # H-alpha
    'OIII': {'central_wavelength_nm': 501, 'color_rgb': (0.0, 1.0, 0.4)}, # [OIII]
}

def get_filter_color(filter_name: str) -> Tuple[float, float, float]:
    """
    Get RGB color for astronomical filter
    
    Args:
        filter_name: Name of the filter
    
    Returns:
        RGB color tuple
    """
    return COMMON_FILTERS.get(filter_name, (1.0, 1.0, 1.0))['color_rgb']
