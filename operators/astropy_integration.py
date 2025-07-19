"""
Comprehensive Astropy Integration
================================

Complete integration of Astropy core modules for professional astronomical analysis.
Covers WCS, coordinates, units, constants, time, cosmology, and IO operations.
"""

import bpy
import numpy as np
from pathlib import Path
from bpy.props import StringProperty, FloatProperty, BoolProperty, EnumProperty, IntProperty, FloatVectorProperty
from bpy_extras.io_utils import ImportHelper

# Core Astropy imports - following standards (no try/except)
from astropy import units as u
from astropy import constants as const
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS, FK4, FK5, Galactic
from astropy.time import Time
from astropy.cosmology import FlatLambdaCDM, Planck18, WMAP9, z_at_value
from astropy.wcs import WCS
from astropy.io import fits
from astropy.table import Table

from ..utilities.astronomical_data import get_stellar_color_from_temperature


class ASTRO_OT_coordinate_converter(bpy.types.Operator):
    """Convert between different astronomical coordinate systems"""
    
    bl_idname = "astronomical.coordinate_converter"
    bl_label = "Coordinate Converter"
    bl_description = "Convert between astronomical coordinate systems using Astropy"
    bl_options = {'REGISTER', 'UNDO'}
    
    # Input coordinates
    input_ra: FloatProperty(
        name="RA (degrees)",
        description="Right Ascension in degrees",
        default=266.405,  # Galactic center
        min=0.0,
        max=360.0
    )
    
    input_dec: FloatProperty(
        name="Dec (degrees)", 
        description="Declination in degrees",
        default=-29.006,  # Galactic center
        min=-90.0,
        max=90.0
    )
    
    input_frame: EnumProperty(
        name="Input Frame",
        description="Input coordinate frame",
        items=[
            ('ICRS', "ICRS", "International Celestial Reference System"),
            ('FK4', "FK4", "Fourth Fundamental Catalogue"),
            ('FK5', "FK5", "Fifth Fundamental Catalogue"),
            ('GALACTIC', "Galactic", "Galactic coordinates"),
            ('ECLIPTIC', "Ecliptic", "Ecliptic coordinates")
        ],
        default='ICRS'
    )
    
    output_frame: EnumProperty(
        name="Output Frame",
        description="Output coordinate frame", 
        items=[
            ('ICRS', "ICRS", "International Celestial Reference System"),
            ('FK4', "FK4",