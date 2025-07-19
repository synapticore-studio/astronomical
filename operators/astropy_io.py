"""
Astropy I/O Operators
====================

Advanced I/O operators for astronomical data formats using Astropy.
Supports FITS, ASCII tables, VOTables, coordinates, and time handling.
"""

import bpy
import numpy as np
from pathlib import Path
from bpy.props import StringProperty, BoolProperty, EnumProperty, FloatProperty
from bpy_extras.io_utils import ImportHelper, ExportHelper

# Astropy imports
import astropy.units as u
import astropy.coordinates as coord
import astropy.time as time
import astropy.io.fits as fits
import astropy.io.ascii as ascii_io
import astropy.io.votable as votable
import astropy.table as table
from astropy.coordinates import SkyCoord
from astropy.time import Time

from ..utilities.astropy_integration import (
    AstronomicalConstants, CoordinateSystem, AstronomicalTime,
    WorldCoordinateSystem, AstropyIO, UnitConversions, TableOperations
)


class ASTRO_OT_import_table(bpy.types.Operator, ImportHelper):
    """Import astronomical table data (FITS, ASCII, VOTable)"""
    
    bl_idname = "astronomical.import_table"
    bl_label = "Import Astronomical Table"
    bl_description = "Import astronomical catalogs and tables using Astropy"
    bl_options = {'REGISTER', 'UNDO'}
    
    # File filter
    filename_ext = ".fits"
    filter_glob: StringProperty(
        default="*.fits;*.fit;*.csv;*.txt;*.dat;*.xml;*.vot",
        options={'HIDDEN'}
    )
    
    # Table format
    table_format: EnumProperty(
        name="Table Format",
        description="Format of the table file",
        items=[
            ('AUTO', "Auto-detect", "Automatically detect format"),
            ('FITS', "FITS Table", "FITS binary or ASCII table"),
            ('ASCII', "ASCII Table", "Generic ASCII table"),
            ('CSV', "CSV", "Comma-separated values"),
            ('VOTABLE', "VOTable", "Virtual Observatory Table"),
            ('IPAC', "IPAC Table", "IPAC ASCII table format")
        ],
        default='AUTO'
    )
    
    # Processing options
    create_skycoord: BoolProperty(
        name="Create Sky Coordinates",
        description="Create SkyCoord objects from RA/Dec columns",
        default=True
    )
    
    ra_column: StringProperty(
        name="RA Column",
        description="Name of RA column (auto-detected if empty)",
        default=""
    )
    
    dec_column: StringProperty(
        name="Dec Column", 
        description="Name of Dec column (auto-detected if empty)",
        default=""
    )
    
    coordinate_unit: EnumProperty(
        name="Coordinate Unit",
        description="Unit for coordinates",
        items=[
            ('deg', "Degrees", "Decimal degrees"),
            ('hour', "Hours", "RA in hours, Dec in degrees"),
            ('rad', "Radians", "Radians")
        ],
        default='deg'
    )
    
    # Visualization options
    create_point_cloud: BoolProperty(
        name="Create Point Cloud",
        description="Create 3D point cloud visualization",
        default=False
    )
    
    distance_column: StringProperty(
        name="Distance Column",
        description="Column for distance/redshift (for 3D visualization)",
        default=""
    )
    
    magnitude_column: StringProperty(
        name="Magnitude Column",
        description="Column for magnitude/brightness",
        default=""
    )

    def execute(self, context):
        try:
            # Import table
            table_data = self.import_table_data()
            
            # Process coordinates if requested
            if self.create_skycoord:
                table_data = self.add_coordinate_columns(table_data)
            
            # Create visualizations
            results = self.create_visualizations(table_data, context)
            
            # Store table data for further analysis
            self.store_table_metadata(table_data, context)
            
            self.report({'INFO'}, f"Imported table with {len(table_data)} rows")
            return {'FINISHED'}
            
        except Exception as e:
            self.report({'ERROR'}, f"Table import failed: {str(e)}")
            return {'CANCELLED'}
    
    def import_table_data(self):
        """Import table using appropriate Astropy reader"""
        
        filepath = Path(self.filepath)
        
        if self.table_format == 'AUTO':
            # Auto-detect format based on extension
            ext = filepath.suffix.lower()
            if ext in ['.fits', '.fit']:
                format_type = 'fits'
            elif ext in ['.xml', '.vot']:
                format_type = 'votable'
            elif ext == '.csv':
                format_type = 'csv'
            else:
                format_type = 'ascii'
        else:
            format_type = self.table_format.lower()
        
        # Import with appropriate reader
        if format_type == 'fits':
            table_data = AstropyIO.read_fits_table(filepath)
        elif format_type == 'votable':
            table_data = AstropyIO.read_votable(filepath)
        elif format_type == 'csv':
            table_data = AstropyIO.read_ascii_table(filepath, format='csv')
        else:
            table_data = AstropyIO.read_ascii_table(filepath, format='basic')
        
        return table_data
    
    def add_coordinate_columns(self, table_data):
        """Add SkyCoord columns to table"""
        
        # Auto-detect RA/Dec columns if not specified
        ra_col = self.ra_column
        dec_col = self.dec_column
        
        if not ra_col:
            # Common RA column names
            ra_candidates = ['RA', 'ra', 'RAJ2000', 'RA_J2000', 'alpha', 'ALPHA_J2000']
            for candidate in ra_candidates:
                if candidate in table_data.colnames:
                    ra_col = candidate
                    break
        
        if not dec_col:
            # Common Dec column names
            dec_candidates = ['DEC', 'dec', 'DEJ2000', 'DE_J2000', 'delta', 'DELTA_J2000']
            for candidate in dec_candidates:
                if candidate in table_data.colnames:
                    dec_col = candidate
                    break
        
        if ra_col and dec_col:
            # Create SkyCoord column
            table_data = TableOperations.add_skycoord_columns(
                table_data, ra_col, dec_col, 'skycoord'
            )
            
            # Add galactic coordinates
            galactic_coords = table_data['skycoord'].galactic
            table_data['galactic_l'] = galactic_coords.l
            table_data['galactic_b'] = galactic_coords.b
            
            self.report({'INFO'}, f"Created coordinates from {ra_col}, {dec_col}")
        else:
            self.report({'WARNING'}, "Could not find RA/Dec columns for coordinates")
        
        return table_data
    
    def create_visualizations(self, table_data, context):
        """Create Blender visualizations from table data"""
        
        results = {}
        
        if self.create_point_cloud and 'skycoord' in table_data.colnames:
            point_cloud = self.create_astronomical_point_cloud(table_data, context)
            results['point_cloud'] = point_cloud
        
        return results
    
    def create_astronomical_point_cloud(self, table_data, context):
        """Create 3D point cloud from astronomical catalog"""
        
        # Get coordinates
        coords = table_data['skycoord']
        ra_deg = coords.ra.deg
        dec_deg = coords.dec.deg
        
        # Get distances if available
        if self.distance_column and self.distance_column in table_data.colnames:
            distances = table_data[self.distance_column]
            # Convert to physical coordinates
            x = distances * np.cos(np.radians(dec_deg)) * np.cos(np.radians(ra_deg))
            y = distances * np.cos(np.radians(dec_deg)) * np.sin(np.radians(ra_deg))
            z = distances * np.sin(np.radians(dec_deg))
        else:
            # Project onto unit sphere
            x = np.cos(np.radians(dec_deg)) * np.cos(np.radians(ra_deg))
            y = np.cos(np.radians(dec_deg)) * np.sin(np.radians(ra_deg))
            z = np.sin(np.radians(dec_deg))
        
        # Scale for Blender
        scale_factor = 10.0
        vertices = [(x[i] * scale_factor, y[i] * scale_factor, z[i] * scale_factor) 
                   for i in range(len(x))]
        
        # Create mesh
        mesh = bpy.data.meshes.new("AstronomicalCatalog")
        mesh.from_pydata(vertices, [], [])
        
        # Create object
        obj = bpy.data.objects.new("AstronomicalCatalog", mesh)
        context.collection.objects.link(obj)
        
        # Create material based on magnitude if available
        if self.magnitude_column and self.magnitude_column in table_data.colnames:
            self.create_magnitude_material(obj, table_data[self.magnitude_column])
        
        return obj
    
    def create_magnitude_material(self, obj, magnitudes):
        """Create material with magnitude-based colors"""
        
        # Create vertex colors based on magnitude
        mesh = obj.data
        mesh.vertex_colors.new()
        color_layer = mesh.vertex_colors.active
        
        # Normalize magnitudes to colors
        mag_min, mag_max = np.min(magnitudes), np.max(magnitudes)
        
        for i, poly in enumerate(mesh.polygons):
            for loop_idx in poly.loop_indices:
                vert_idx = mesh.loops[loop_idx].vertex_index
                if vert_idx < len(magnitudes):
                    # Bright objects (low magnitude) = white, faint = red
                    norm_mag = (magnitudes[vert_idx] - mag_min) / (mag_max - mag_min)
                    color = (1.0 - norm_mag, 0.5, norm_mag, 1.0)
                    color_layer.data[loop_idx].color = color
    
    def store_table_metadata(self, table_data, context):
        """Store table metadata in scene"""
        
        # Store basic info as scene properties
        context.scene["astropy_table_imported"] = True
        context.scene["astropy_table_rows"] = len(table_data)
        context.scene["astropy_table_columns"] = len(table_data.colnames)
        context.scene["astropy_table_source"] = str(self.filepath)


class ASTRO_OT_coordinate_transform(bpy.types.Operator):
    """Transform coordinates between different systems"""
    
    bl_idname = "astronomical.coordinate_transform" 
    bl_label = "Coordinate Transform"
    bl_description = "Transform coordinates between astronomical coordinate systems"
    bl_options = {'REGISTER', 'UNDO'}
    
    # Input coordinates
    input_ra: FloatProperty(
        name="RA / Longitude",
        description="Right Ascension or Longitude",
        default=0.0
    )
    
    input_dec: FloatProperty(
        name="Dec / Latitude", 
        description="Declination or Latitude",
        default=0.0
    )
    
    input_frame: EnumProperty(
        name="Input Frame",
        description="Input coordinate frame",
        items=[
            ('icrs', "ICRS", "International Celestial Reference System"),
            ('fk5', "FK5", "FK5 (J2000.0)"),
            ('galactic', "Galactic", "Galactic coordinates"),
            ('ecliptic', "Ecliptic", "Ecliptic coordinates"),
            ('supergalactic', "Supergalactic", "Supergalactic coordinates")
        ],
        default='icrs'
    )
    
    output_frame: EnumProperty(
        name="Output Frame",
        description="Output coordinate frame", 
        items=[
            ('icrs', "ICRS", "International Celestial Reference System"),
            ('fk5', "FK5", "FK5 (J2000.0)"),
            ('galactic', "Galactic", "Galactic coordinates"),
            ('ecliptic', "Ecliptic", "Ecliptic coordinates"),
            ('supergalactic', "Supergalactic", "Supergalactic coordinates")
        ],
        default='galactic'
    )
    
    coordinate_unit: EnumProperty(
        name="Input Unit",
        description="Unit for input coordinates",
        items=[
            ('deg', "Degrees", "Decimal degrees"),
            ('hour', "Hours", "RA in hours, Dec in degrees")
        ],
        default='deg'
    )

    def execute(self, context):
        try:
            # Create input coordinates
            if self.coordinate_unit == 'hour':
                coord_obj = CoordinateSystem.create_skycoord(
                    self.input_ra, self.input_dec, frame=self.input_frame, unit='hour'
                )
            else:
                coord_obj = CoordinateSystem.create_skycoord(
                    self.input_ra, self.input_dec, frame=self.input_frame, unit='deg'
                )
            
            # Transform to output frame
            transformed = CoordinateSystem.transform_coordinates(coord_obj, self.output_frame)
            
            # Format results
            if self.output_frame == 'galactic':
                result_text = f"l = {transformed.l.deg:.6f}°, b = {transformed.b.deg:.6f}°"
            elif self.output_frame == 'ecliptic':
                result_text = f"lon = {transformed.lon.deg:.6f}°, lat = {transformed.lat.deg:.6f}°"
            else:
                result_text = f"RA = {transformed.ra.deg:.6f}°, Dec = {transformed.dec.deg:.6f}°"
                result_text += f"\nRA = {transformed.ra.to_string(unit=u.hour)}, Dec = {transformed.dec.to_string(unit=u.deg)}"
            
            # Store results
            context.scene["coordinate_transform_result"] = result_text
            
            self.report({'INFO'}, f"Transformed: {result_text}")
            return {'FINISHED'}
            
        except Exception as e:
            self.report({'ERROR'}, f"Coordinate transformation failed: {str(e)}")
            return {'CANCELLED'}

    def draw(self, context):
        layout = self.layout
        
        layout.prop(self, "input_frame")
        layout.prop(self, "coordinate_unit")
        
        row = layout.row()
        row.prop(self, "input_ra")
        row.prop(self, "input_dec")
        
        layout.prop(self, "output_frame")
        
        # Show result if available
        if "coordinate_transform_result" in context.scene:
            box = layout.box()
            box.label(text="Result:")
            for line in context.scene["coordinate_transform_result"].split('\n'):
                box.label(text=line)


class ASTRO_OT_time_converter(bpy.types.Operator):
    """Convert between different time systems"""
    
    bl_idname = "astronomical.time_converter"
    bl_label = "Time Converter"
    bl_description = "Convert between astronomical time systems"
    bl_options = {'REGISTER', 'UNDO'}
    
    # Input time
    input_time: StringProperty(
        name="Input Time",
        description="Input time (ISO format, MJD, or JD)",
        default="2024-01-01T00:00:00"
    )
    
    input_format: EnumProperty(
        name="Input Format",
        description="Format of input time",
        items=[
            ('iso', "ISO", "ISO 8601 format (YYYY-MM-DDTHH:MM:SS)"),
            ('mjd', "MJD", "Modified Julian Date"),
            ('jd', "JD", "Julian Date"),
            ('fits', "FITS", "FITS time format")
        ],
        default='iso'
    )
    
    input_scale: EnumProperty(
        name="Input Scale",
        description="Time scale of input",
        items=[
            ('utc', "UTC", "Coordinated Universal Time"),
            ('tai', "TAI", "International Atomic Time"),
            ('tt', "TT", "Terrestrial Time"),
            ('tdb', "TDB", "Barycentric Dynamical Time")
        ],
        default='utc'
    )
    
    output_format: EnumProperty(
        name="Output Format",
        description="Desired output format",
        items=[
            ('iso', "ISO", "ISO 8601 format"),
            ('mjd', "MJD", "Modified Julian Date"),
            ('jd', "JD", "Julian Date"),
            ('fits', "FITS", "FITS time format")
        ],
        default='mjd'
    )
    
    output_scale: EnumProperty(
        name="Output Scale", 
        description="Desired output time scale",
        items=[
            ('utc', "UTC", "Coordinated Universal Time"),
            ('tai', "TAI", "International Atomic Time"), 
            ('tt', "TT", "Terrestrial Time"),
            ('tdb', "TDB", "Barycentric Dynamical Time")
        ],
        default='utc'
    )

    def execute(self, context):
        try:
            # Create Time object
            time_obj = AstronomicalTime.create_time(
                self.input_time, format=self.input_format, scale=self.input_scale
            )
            
            # Convert to output scale
            if self.output_scale != self.input_scale:
                time_obj = getattr(time_obj, self.output_scale)
            
            # Format output
            if self.output_format == 'iso':
                result = AstronomicalTime.time_to_string(time_obj, 'iso')
            elif self.output_format == 'mjd':
                result = f"{time_obj.mjd:.6f}"
            elif self.output_format == 'jd':
                result = f"{time_obj.jd:.6f}"
            elif self.output_format == 'fits':
                result = AstronomicalTime.time_to_string(time_obj, 'fits')
            
            # Store result
            context.scene["time_conversion_result"] = result
            
            self.report({'INFO'}, f"Converted time: {result}")
            return {'FINISHED'}
            
        except Exception as e:
            self.report({'ERROR'}, f"Time conversion failed: {str(e)}")
            return {'CANCELLED'}

    def draw(self, context):
        layout = self.layout
        
        layout.prop(self, "input_time")
        row = layout.row()
        row.prop(self, "input_format")
        row.prop(self, "input_scale")
        
        layout.separator()
        
        row = layout.row()
        row.prop(self, "output_format")
        row.prop(self, "output_scale")
        
        # Show result if available
        if "time_conversion_result" in context.scene:
            box = layout.box()
            box.label(text=f"Result: {context.scene['time_conversion_result']}")


class ASTRO_OT_export_table(bpy.types.Operator, ExportHelper):
    """Export astronomical table data"""
    
    bl_idname = "astronomical.export_table"
    bl_label = "Export Astronomical Table"
    bl_description = "Export astronomical data as table"
    bl_options = {'REGISTER', 'UNDO'}
    
    filename_ext = ".fits"
    filter_glob: StringProperty(
        default="*.fits;*.csv;*.txt;*.xml",
        options={'HIDDEN'}
    )
    
    export_format: EnumProperty(
        name="Export Format",
        description="Format for exported table",
        items=[
            ('FITS', "FITS Table", "FITS binary table"),
            ('CSV', "CSV", "Comma-separated values"),
            ('ASCII', "ASCII", "ASCII table"),
            ('VOTABLE', "VOTable", "Virtual Observatory Table")
        ],
        default='FITS'
    )

    @classmethod
    def poll(cls, context):
        return context.scene.get("astropy_table_imported", False)

    def execute(self, context):
        try:
            # This would export data from stored table
            # For now, create a simple example
            self.report({'INFO'}, "Table export functionality - implementation depends on stored data structure")
            return {'FINISHED'}
            
        except Exception as e:
            self.report({'ERROR'}, f"Table export failed: {str(e)}")
            return {'CANCELLED'}


# Panel for Astropy I/O tools
class ASTRO_PT_astropy_io_panel(bpy.types.Panel):
    """Astropy I/O tools panel"""
    
    bl_label = "Astropy I/O"
    bl_idname = "ASTRO_PT_astropy_io"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Astronomical"
    bl_parent_id = "ASTRO_PT_astrophot"
    
    def draw(self, context):
        layout = self.layout
        
        # Table I/O
        box = layout.box()
        box.label(text="Table I/O:", icon='SPREADSHEET')
        box.operator("astronomical.import_table", text="Import Table", icon='IMPORT')
        if context.scene.get("astropy_table_imported", False):
            box.operator("astronomical.export_table", text="Export Table", icon='EXPORT')
            box.label(text=f"Rows: {context.scene.get('astropy_table_rows', 0)}")
        
        # Coordinate tools
        box = layout.box()
        box.label(text="Coordinates:", icon='ORIENTATION_GLOBAL')
        box.operator("astronomical.coordinate_transform", text="Transform Coordinates")
        
        # Time tools
        box = layout.box()
        box.label(text="Time Systems:", icon='TIME')
        box.operator("astronomical.time_converter", text="Convert Time")
        
        # Constants info
        box = layout.box()
        box.label(text="Constants Available:", icon='INFO')
        constants = AstronomicalConstants.get_constant_dict()
        box.label(text=f"c = {constants['c']:.2e}")
        box.label(text=f"G = {constants['G']:.2e}")
        box.label(text=f"M_sun = {constants['M_sun']:.2e}")


# Registration
classes = [
    ASTRO_OT_import_table,
    ASTRO_OT_coordinate_transform,
    ASTRO_OT_time_converter,
    ASTRO_OT_export_table,
    ASTRO_PT_astropy_io_panel,
]

def register():
    """Register Astropy I/O operators"""
    for cls in classes:
        bpy.utils.register_class(cls)

def unregister():
    """Unregister Astropy I/O operators"""
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)

if __name__ == "__main__":
    register()
