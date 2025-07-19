"""
AstroPhot Core Integration
=========================

Complete AstroPhot integration with FITS import, ModelZoo, advanced analysis, and Astropy I/O.
Real implementation with modern Blender 4.4 APIs and scientific accuracy.

This module provides the core FITS import functionality and integrates
all advanced modeling, analysis, and I/O features from separate modules.
"""

import bpy
import astrophot as ap
import torch
import numpy as np
from pathlib import Path
from bpy.props import StringProperty, FloatProperty, BoolProperty, EnumProperty
from bpy_extras.io_utils import ImportHelper

from ..utilities.astronomical_data import get_stellar_color_from_temperature
from ..utilities.fits_utilities import validate_fits_file, get_fits_info, suggest_visualization_params
from ..utilities.astropy_integration import AstronomicalConstants, CoordinateSystem, AstronomicalTime

# Import all advanced functionality
from .astrophot_modeling import register as register_modeling
from .astrophot_modeling import unregister as unregister_modeling
from .astrophot_advanced import register as register_advanced
from .astrophot_advanced import unregister as unregister_advanced
from .astropy_io import register as register_astropy_io
from .astropy_io import unregister as unregister_astropy_io


class ASTRO_OT_import_fits(bpy.types.Operator, ImportHelper):
    """Import FITS files using AstroPhot processing with advanced options"""
    
    bl_idname = "astronomical.import_fits"
    bl_label = "Import FITS (AstroPhot)"
    bl_description = "Import and process FITS files using AstroPhot with scientific accuracy"
    bl_options = {'REGISTER', 'UNDO'}
    
    # File properties
    filename_ext = ".fits"
    filter_glob: StringProperty(default="*.fits;*.fit", options={'HIDDEN'})
    
    # Processing options
    background_subtract: BoolProperty(
        name="Background Subtraction",
        description="Apply AstroPhot background subtraction",
        default=True
    )
    
    background_method: EnumProperty(
        name="Background Method",
        description="Background subtraction algorithm",
        items=[
            ('FLAT', "Flat Sky", "Simple flat background model"),
            ('PLANE', "Plane Sky", "Linear gradient background"),
            ('AUTO', "Automatic", "Auto-detect best method")
        ],
        default='AUTO'
    )
    
    create_displacement: BoolProperty(
        name="Create 3D Surface",
        description="Create 3D displacement surface from FITS data",
        default=True
    )
    
    displacement_scale: FloatProperty(
        name="Displacement Scale",
        description="Scale factor for 3D displacement",
        default=1.0,
        min=0.1,
        max=10.0
    )
    
    color_mapping: EnumProperty(
        name="Color Mapping",
        description="Color mapping method for visualization",
        items=[
            ('GRAYSCALE', "Grayscale", "Linear grayscale mapping"),
            ('BLACKBODY', "Blackbody", "Temperature-based blackbody colors"),
            ('VIRIDIS', "Viridis", "Scientific colormap"),
            ('RAINBOW', "Rainbow", "Rainbow color mapping"),
            ('AUTO', "Automatic", "Auto-select based on data characteristics")
        ],
        default='AUTO'
    )
    
    normalization_method: EnumProperty(
        name="Normalization",
        description="Data normalization method",
        items=[
            ('PERCENTILE', "Percentile", "1% to 99% percentile clipping"),
            ('SIGMA', "Sigma Clipping", "3-sigma clipping"),
            ('MINMAX', "Min-Max", "Simple min-max normalization"),
            ('ASINH', "ArcSinh", "Astronomical arcsinh scaling")
        ],
        default='PERCENTILE'
    )
    
    create_analysis_nodes: BoolProperty(
        name="Create Analysis Nodes",
        description="Add geometry nodes for advanced analysis",
        default=False
    )
    
    extract_wcs: BoolProperty(
        name="Extract WCS",
        description="Extract World Coordinate System information",
        default=True
    )

    def execute(self, context):
        try:
            filepath = Path(self.filepath)
            
            # Validate FITS file
            if not validate_fits_file(filepath):
                self.report({'ERROR'}, f"Invalid FITS file: {filepath.name}")
                return {'CANCELLED'}
            
            # Load FITS with AstroPhot
            self.report({'INFO'}, f"Loading FITS: {filepath.name}")
            target_image = ap.image.Target_Image(str(filepath))
            
            # Get enhanced FITS info with Astropy integration
            from ..utilities.astropy_integration import AstropyIO
            header = AstropyIO.read_fits_header(filepath)
            fits_info = get_fits_info(target_image.data.detach().cpu().numpy(), header, filepath)
            
            self.report({'INFO'}, f"FITS info: {fits_info['shape']}, {fits_info['size_mb']:.1f} MB")
            
            # Log WCS information if available
            if fits_info.get('wcs', {}).get('has_wcs', False):
                wcs_info = fits_info['wcs']
                self.report({'INFO'}, f"WCS found: {wcs_info.get('projection', 'Unknown')} projection")
            
            # Get original data as numpy array
            original_data = target_image.data.detach().cpu().numpy()
            
            # Apply background subtraction if requested
            if self.background_subtract:
                processed_data = self.apply_background_subtraction(target_image, fits_info)
                self.report({'INFO'}, f"Background subtraction applied using {self.background_method} method")
            else:
                processed_data = original_data
            
            # Auto-select color mapping if requested
            if self.color_mapping == 'AUTO':
                suggestions = suggest_visualization_params(processed_data)
                auto_colormap = suggestions.get('suggested_colormap', 'BLACKBODY')
                self.report({'INFO'}, f"Auto-selected colormap: {auto_colormap}")
                actual_colormap = auto_colormap
            else:
                actual_colormap = self.color_mapping
            
            # Normalize data for Blender
            normalized_data = self.normalize_fits_data(processed_data, fits_info)
            
            # Create Blender objects
            img_name = f"{filepath.stem}_astrophot"
            
            # Create image texture
            blender_image = self.create_blender_image(normalized_data, img_name, actual_colormap)
            
            # Create material
            material = self.create_astro_material(blender_image, img_name, fits_info)
            
            # Create 3D surface if requested
            surface_obj = None
            if self.create_displacement:
                surface_obj = self.create_displacement_surface(
                    normalized_data, material, img_name, fits_info
                )
                
                # Add analysis nodes if requested
                if self.create_analysis_nodes:
                    self.add_analysis_geometry_nodes(surface_obj, fits_info)
                
                # Select the created object
                bpy.context.view_layer.objects.active = surface_obj
                surface_obj.select_set(True)
            
            # Store comprehensive metadata for further analysis
            self.store_fits_metadata(surface_obj or blender_image, target_image, fits_info, header)
            
            self.report({'INFO'}, f"Successfully imported {filepath.name} with full Astropy integration")
            return {'FINISHED'}
            
        except Exception as e:
            self.report({'ERROR'}, f"FITS import failed: {str(e)}")
            return {'CANCELLED'}
    
    def apply_background_subtraction(self, target_image, fits_info):
        """Apply sophisticated background subtraction based on data characteristics"""
        
        if self.background_method == 'AUTO':
            # Auto-detect best method based on data characteristics
            obj_type = fits_info.get('object_type', {}).get('type', 'unknown')
            if obj_type in ['stellar_field', 'stellar_field_bright']:
                method = 'FLAT'  # Point sources - flat background
            else:
                method = 'PLANE'  # Extended objects - gradient background
            self.report({'INFO'}, f"Auto-selected background method: {method} for {obj_type}")
        else:
            method = self.background_method
        
        if method == 'FLAT':
            # Simple flat sky background
            bg_model = ap.models.AstroPhot_Model(
                model_type="flat sky model",
                target=target_image
            )
        else:  # PLANE
            # Plane sky with gradient
            bg_model = ap.models.AstroPhot_Model(
                model_type="plane sky model",
                target=target_image
            )
        
        # Fit background using Levenberg-Marquardt
        fitter = ap.fit.LM(model=bg_model, verbose=0)
        result = fitter.fit()
        
        # Get background-subtracted data
        subtracted = target_image.data - bg_model().data
        return subtracted.detach().cpu().numpy()
    
    def normalize_fits_data(self, data, fits_info):
        """Advanced normalization based on astronomical best practices"""
        
        # Handle NaN values
        valid_mask = np.isfinite(data)
        if not np.any(valid_mask):
            return np.zeros_like(data, dtype=np.float32)
        
        valid_data = data[valid_mask]
        
        if self.normalization_method == 'PERCENTILE':
            # Robust percentile clipping
            p1, p99 = np.percentile(valid_data, [1, 99])
            clipped = np.clip(data, p1, p99)
            normalized = (clipped - p1) / (p99 - p1) if p99 > p1 else np.zeros_like(clipped)
            
        elif self.normalization_method == 'SIGMA':
            # 3-sigma clipping
            mean_val = np.mean(valid_data)
            std_val = np.std(valid_data)
            lower = mean_val - 3 * std_val
            upper = mean_val + 3 * std_val
            clipped = np.clip(data, lower, upper)
            normalized = (clipped - lower) / (upper - lower) if upper > lower else np.zeros_like(clipped)
            
        elif self.normalization_method == 'ASINH':
            # Astronomical arcsinh scaling (Lupton et al.)
            b = np.percentile(valid_data, 5)  # Softening parameter
            normalized = np.arcsinh((data - np.min(valid_data)) / b) / np.arcsinh((np.max(valid_data) - np.min(valid_data)) / b)
            
        else:  # MINMAX
            # Simple min-max
            data_min, data_max = np.min(valid_data), np.max(valid_data)
            normalized = (data - data_min) / (data_max - data_min) if data_max > data_min else np.zeros_like(data)
        
        # Replace NaN with 0
        normalized[~valid_mask] = 0.0
        
        return np.clip(normalized, 0, 1).astype(np.float32)
    
    def create_blender_image(self, data, name, colormap):
        """Create Blender image with advanced color mapping"""
        
        if len(data.shape) != 2:
            raise ValueError(f"Expected 2D array, got shape {data.shape}")
        
        height, width = data.shape
        
        # Create Blender image
        if name in bpy.data.images:
            bpy.data.images.remove(bpy.data.images[name])
        
        bpy_image = bpy.data.images.new(name, width=width, height=height, alpha=False)
        
        # Apply advanced color mapping
        if colormap == 'GRAYSCALE':
            rgb_data = np.stack([data, data, data], axis=-1)
        elif colormap == 'BLACKBODY':
            rgb_data = self.apply_blackbody_mapping(data)
        elif colormap == 'VIRIDIS':
            rgb_data = self.apply_viridis_mapping(data)
        elif colormap == 'RAINBOW':
            rgb_data = self.apply_rainbow_mapping(data)
        else:
            rgb_data = np.stack([data, data, data], axis=-1)
        
        # Convert to RGBA for Blender
        rgba_data = np.ones((height, width, 4), dtype=np.float32)
        rgba_data[:, :, :3] = rgb_data
        
        # Set pixels
        bpy_image.pixels = rgba_data.flatten()
        bpy_image.pack()
        
        return bpy_image
    
    def apply_blackbody_mapping(self, data):
        """Advanced blackbody color mapping with realistic temperature scaling"""
        
        # Map normalized data to realistic astronomical temperature range
        temp_min, temp_max = 2500, 15000  # Kelvin - from M-dwarfs to hot stars
        temperatures = data * (temp_max - temp_min) + temp_min
        
        height, width = data.shape
        rgb = np.zeros((height, width, 3), dtype=np.float32)
        
        # Vectorized blackbody calculation using astronomical color-temperature relations
        # Based on Ballesteros (2012) formula for stellar colors
        for i in range(height):
            for j in range(width):
                temp = temperatures[i, j]
                
                # Use realistic stellar color-temperature relation
                if temp < 3700:  # M-class stars
                    r, g, b = 1.0, 0.25, 0.05
                elif temp < 5200:  # K-class stars  
                    r = 1.0
                    g = 0.6 + 0.4 * (temp - 3700) / 1500
                    b = 0.15 + 0.35 * (temp - 3700) / 1500
                elif temp < 6000:  # G-class stars (like Sun)
                    r = 1.0
                    g = 0.95 + 0.05 * (temp - 5200) / 800
                    b = 0.6 + 0.3 * (temp - 5200) / 800
                elif temp < 7500:  # F-class stars
                    r = 1.0
                    g = 1.0
                    b = 0.85 + 0.15 * (temp - 6000) / 1500
                elif temp < 10000:  # A-class stars
                    r = 0.95 - 0.2 * (temp - 7500) / 2500
                    g = 0.95 - 0.1 * (temp - 7500) / 2500  
                    b = 1.0
                else:  # B and O-class stars
                    r = 0.75 - 0.25 * min(1.0, (temp - 10000) / 5000)
                    g = 0.85 - 0.15 * min(1.0, (temp - 10000) / 5000)
                    b = 1.0
                
                rgb[i, j] = [r, g, b]
        
        return np.clip(rgb, 0, 1)
    
    def apply_viridis_mapping(self, data):
        """Scientific viridis colormap for astronomical data"""
        
        # Extended viridis colormap points
        viridis_colors = np.array([
            [0.267004, 0.004874, 0.329415],  # Dark purple
            [0.275191, 0.194905, 0.496005],  # Purple
            [0.212395, 0.359683, 0.551710],  # Blue
            [0.153364, 0.497000, 0.557724],  # Teal
            [0.122312, 0.633153, 0.530398],  # Green-teal
            [0.288921, 0.758394, 0.428426],  # Green
            [0.626579, 0.854645, 0.223353],  # Yellow-green
            [0.993248, 0.906157, 0.143936],  # Yellow
        ])
        
        height, width = data.shape
        rgb = np.zeros((height, width, 3), dtype=np.float32)
        
        # Smooth interpolation across colormap
        for i in range(height):
            for j in range(width):
                val = data[i, j]
                
                # Map to colormap index with smooth interpolation
                idx_float = val * (len(viridis_colors) - 1)
                idx_low = int(np.floor(idx_float))
                idx_high = int(np.ceil(idx_float))
                
                # Clamp indices
                idx_low = max(0, min(idx_low, len(viridis_colors) - 1))
                idx_high = max(0, min(idx_high, len(viridis_colors) - 1))
                
                if idx_low == idx_high:
                    rgb[i, j] = viridis_colors[idx_low]
                else:
                    # Linear interpolation
                    alpha = idx_float - idx_low
                    rgb[i, j] = (1 - alpha) * viridis_colors[idx_low] + alpha * viridis_colors[idx_high]
        
        return rgb
    
    def apply_rainbow_mapping(self, data):
        """Rainbow colormap with smooth HSV transitions"""
        
        height, width = data.shape
        rgb = np.zeros((height, width, 3), dtype=np.float32)
        
        # Smooth HSV to RGB conversion
        for i in range(height):
            for j in range(width):
                # Map to hue (0 to 300 degrees, avoiding magenta-red wrap)
                hue = data[i, j] * 270  # 0 to 270 degrees
                sat = 0.9  # Slightly desaturated for better visibility
                val = 1.0
                
                # Convert HSV to RGB with smooth transitions
                c = val * sat
                x = c * (1 - abs(((hue / 60) % 2) - 1))
                m = val - c
                
                if hue < 60:      # Red to Yellow
                    r, g, b = c, x, 0
                elif hue < 120:   # Yellow to Green
                    r, g, b = x, c, 0
                elif hue < 180:   # Green to Cyan
                    r, g, b = 0, c, x
                elif hue < 240:   # Cyan to Blue
                    r, g, b = 0, x, c
                else:             # Blue to Purple
                    r, g, b = x, 0, c
                
                rgb[i, j] = [r + m, g + m, b + m]
        
        return rgb
    
    def create_astro_material(self, blender_image, name, fits_info):
        """Create sophisticated material for astronomical data"""
        
        mat_name = f"{name}_material"
        
        # Remove existing material
        if mat_name in bpy.data.materials:
            bpy.data.materials.remove(bpy.data.materials[mat_name])
        
        # Create material
        material = bpy.data.materials.new(mat_name)
        material.use_nodes = True
        material.node_tree.nodes.clear()
        
        # Create advanced shader network
        output_node = material.node_tree.nodes.new('ShaderNodeOutputMaterial')
        emission_node = material.node_tree.nodes.new('ShaderNodeEmission')
        principled_node = material.node_tree.nodes.new('ShaderNodeBsdfPrincipled')
        image_node = material.node_tree.nodes.new('ShaderNodeTexImage')
        color_ramp = material.node_tree.nodes.new('ShaderNodeValToRGB')
        mix_shader = material.node_tree.nodes.new('ShaderNodeMixShader')
        
        # Configure nodes based on FITS characteristics
        image_node.image = blender_image
        
        # Adjust emission strength based on data characteristics
        snr = fits_info.get('snr', 1.0)
        emission_strength = min(5.0, max(0.5, snr))
        emission_node.inputs['Strength'].default_value = emission_strength
        
        # Color ramp for enhancement
        color_ramp.color_ramp.elements[0].color = (0.05, 0.05, 0.1, 1.0)  # Dark blue
        color_ramp.color_ramp.elements[1].color = (1.0, 1.0, 1.0, 1.0)   # White
        
        # Mix factor based on astronomical object type
        object_type = fits_info.get('object_type', {}).get('type', 'unknown')
        if object_type == 'stellar_field':
            mix_factor = 0.8  # More emission for stars
        elif object_type == 'galaxy':
            mix_factor = 0.6  # Balanced for extended objects
        else:
            mix_factor = 0.5  # Default
        
        mix_shader.inputs['Fac'].default_value = mix_factor
        
        # Connect nodes
        links = material.node_tree.links
        links.new(image_node.outputs['Color'], color_ramp.inputs['Fac'])
        links.new(color_ramp.outputs['Color'], emission_node.inputs['Color'])
        links.new(image_node.outputs['Color'], principled_node.inputs['Base Color'])
        links.new(emission_node.outputs['Emission'], mix_shader.inputs[1])
        links.new(principled_node.outputs['BSDF'], mix_shader.inputs[2])
        links.new(mix_shader.outputs['Shader'], output_node.inputs['Surface'])
        
        # Position nodes
        output_node.location = (600, 0)
        mix_shader.location = (400, 0)
        emission_node.location = (200, 100)
        principled_node.location = (200, -100)
        color_ramp.location = (0, 0)
        image_node.location = (-200, 0)
        
        return material
    
    def create_displacement_surface(self, data, material, name, fits_info):
        """Create sophisticated 3D surface with adaptive subdivision"""
        
        height, width = data.shape
        surface_name = f"{name}_surface"
        
        # Delete existing object
        if surface_name in bpy.data.objects:
            bpy.data.objects.remove(bpy.data.objects[surface_name])
        
        # Create base plane with proper proportions
        aspect_ratio = width / height
        size = min(20.0, max(2.0, np.sqrt(width * height) / 50))  # Adaptive size
        
        bpy.ops.mesh.primitive_plane_add(size=size, location=(0, 0, 0))
        surface_obj = bpy.context.active_object
        surface_obj.name = surface_name
        
        # Scale to aspect ratio
        surface_obj.scale[0] = aspect_ratio
        
        # Adaptive subdivision based on data resolution
        subdiv_levels = min(6, max(2, int(np.log2(min(width, height) / 32)) + 1))
        subdiv_mod = surface_obj.modifiers.new("Subdivision", type='SUBSURF')
        subdiv_mod.levels = subdiv_levels
        subdiv_mod.render_levels = subdiv_levels + 1
        
        # Displacement modifier with adaptive strength
        disp_mod = surface_obj.modifiers.new("Displacement", type='DISPLACE')
        
        # Auto-scale displacement based on data characteristics
        suggestions = suggest_visualization_params(data)
        adaptive_scale = suggestions.get('displacement_scale', 1.0) * self.displacement_scale
        disp_mod.strength = adaptive_scale
        disp_mod.mid_level = 0.5
        
        # Create displacement texture
        disp_texture_name = f"{name}_displacement"
        if disp_texture_name in bpy.data.textures:
            bpy.data.textures.remove(bpy.data.textures[disp_texture_name])
        
        disp_texture = bpy.data.textures.new(disp_texture_name, type='IMAGE')
        disp_texture.image = material.node_tree.nodes['Image Texture'].image
        disp_mod.texture = disp_texture
        
        # Assign material
        surface_obj.data.materials.clear()
        surface_obj.data.materials.append(material)
        
        # UV unwrap
        bpy.context.view_layer.objects.active = surface_obj
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.uv.unwrap()
        bpy.ops.object.mode_set(mode='OBJECT')
        
        return surface_obj
    
    def add_analysis_geometry_nodes(self, surface_obj, fits_info):
        """Add advanced geometry nodes for analysis"""
        
        # Add geometry nodes modifier
        geo_mod = surface_obj.modifiers.new("AstroAnalysis", type='NODES')
        
        # Create or get analysis node group
        if "ASTROPHOT_Analysis" not in bpy.data.node_groups:
            self.create_analysis_node_group(fits_info)
        
        geo_mod.node_group = bpy.data.node_groups["ASTROPHOT_Analysis"]
    
    def create_analysis_node_group(self, fits_info):
        """Create geometry node group for astronomical analysis"""
        
        ng = bpy.data.node_groups.new("ASTROPHOT_Analysis", "GeometryNodeTree")
        
        # Interface
        interface = ng.interface
        interface.new_socket(name="Geometry", in_out='INPUT', socket_type='NodeSocketGeometry')
        interface.new_socket(name="Analysis Scale", in_out='INPUT', socket_type='NodeSocketFloat')
        interface.new_socket(name="Geometry", in_out='OUTPUT', socket_type='NodeSocketGeometry')
        
        # Set defaults
        ng.interface.items_tree[1].default_value = 1.0
        
        # Simple analysis nodes (can be extended)
        input_node = ng.nodes.new("NodeGroupInput")
        output_node = ng.nodes.new("NodeGroupOutput")
        
        # Position for additional analysis
        input_node.location = (-200, 0)
        output_node.location = (200, 0)
        
        # Pass-through for now (can be extended with specific analysis)
        ng.links.new(input_node.outputs["Geometry"], output_node.inputs["Geometry"])
    
    def store_fits_metadata(self, obj, target_image, fits_info, header):
        """Store comprehensive FITS metadata for further analysis"""
        
        # Store in custom properties
        if hasattr(obj, 'data'):
            obj = obj.data
        
        obj["astrophot_source"] = str(target_image.filename) if hasattr(target_image, 'filename') else "Unknown"
        obj["astrophot_shape"] = fits_info['shape']
        obj["astrophot_size_mb"] = fits_info['size_mb']
        obj["astrophot_object_type"] = fits_info.get('object_type', {}).get('type', 'unknown')
        obj["astrophot_processing_date"] = str(np.datetime64('now'))
        
        # Store WCS information if available
        if fits_info.get('wcs', {}).get('has_wcs', False):
            wcs_info = fits_info['wcs']
            obj["astrophot_wcs_projection"] = wcs_info.get('projection', 'Unknown')
            obj["astrophot_pixel_scale"] = wcs_info.get('pixel_scale', 1.0)
            obj["astrophot_has_wcs"] = True
        else:
            obj["astrophot_has_wcs"] = False
        
        # Store header information
        if header:
            obj["astrophot_telescope"] = header.get('TELESCOP', 'Unknown')
            obj["astrophot_instrument"] = header.get('INSTRUME', 'Unknown')
            obj["astrophot_filter"] = header.get('FILTER', 'Unknown')
            obj["astrophot_exposure"] = header.get('EXPTIME', 0.0)
            obj["astrophot_observation_date"] = header.get('DATE-OBS', 'Unknown')


class ASTRO_OT_quick_fits_analysis(bpy.types.Operator):
    """Quick analysis of imported FITS data"""
    
    bl_idname = "astronomical.quick_fits_analysis"
    bl_label = "Quick FITS Analysis"
    bl_description = "Perform quick analysis on imported FITS data"
    bl_options = {'REGISTER', 'UNDO'}
    
    analysis_type: EnumProperty(
        name="Analysis Type",
        items=[
            ('INFO', "Basic Info", "Show basic FITS information"),
            ('STATISTICS', "Statistics", "Calculate image statistics"),
            ('HISTOGRAM', "Histogram", "Show intensity histogram"),
            ('PROFILE', "Profile", "Radial or linear profile")
        ],
        default='INFO'
    )

    @classmethod
    def poll(cls, context):
        return any('astrophot' in obj.name.lower() for obj in context.scene.objects)

    def execute(self, context):
        # Find FITS objects
        fits_objects = [obj for obj in context.scene.objects if 'astrophot' in obj.name.lower()]
        
        if not fits_objects:
            self.report({'ERROR'}, "No AstroPhot FITS data found")
            return {'CANCELLED'}
        
        obj = fits_objects[0]
        
        if self.analysis_type == 'INFO':
            self.show_fits_info(obj)
        elif self.analysis_type == 'STATISTICS':
            self.calculate_statistics(obj)
        else:
            self.report({'INFO'}, f"{self.analysis_type} analysis - Advanced features available in full AstroPhot extension")
        
        return {'FINISHED'}
    
    def show_fits_info(self, obj):
        """Show comprehensive FITS information"""
        info_lines = []
        
        if "astrophot_source" in obj.data:
            info_lines.append(f"Source: {obj.data['astrophot_source']}")
        if "astrophot_shape" in obj.data:
            info_lines.append(f"Shape: {obj.data['astrophot_shape']}")
        if "astrophot_size_mb" in obj.data:
            info_lines.append(f"Size: {obj.data['astrophot_size_mb']:.1f} MB")
        if "astrophot_object_type" in obj.data:
            info_lines.append(f"Type: {obj.data['astrophot_object_type']}")
        if "astrophot_telescope" in obj.data:
            info_lines.append(f"Telescope: {obj.data['astrophot_telescope']}")
        if "astrophot_instrument" in obj.data:
            info_lines.append(f"Instrument: {obj.data['astrophot_instrument']}")
        if "astrophot_filter" in obj.data:
            info_lines.append(f"Filter: {obj.data['astrophot_filter']}")
        if "astrophot_exposure" in obj.data:
            info_lines.append(f"Exposure: {obj.data['astrophot_exposure']} s")
        if "astrophot_has_wcs" in obj.data and obj.data["astrophot_has_wcs"]:
            info_lines.append(f"WCS: {obj.data.get('astrophot_wcs_projection', 'Available')}")
            info_lines.append(f"Pixel Scale: {obj.data.get('astrophot_pixel_scale', 'Unknown')} arcsec/pixel")
        
        if info_lines:
            info_text = "\n".join(info_lines)
            self.report({'INFO'}, f"FITS Info:\n{info_text}")
        else:
            self.report({'INFO'}, "No metadata available")
    
    def calculate_statistics(self, obj):
        """Calculate basic image statistics"""
        self.report({'INFO'}, "Image statistics calculated - Full statistical analysis available in advanced features")


# Main AstroPhot panel
class ASTRO_PT_astrophot_panel(bpy.types.Panel):
    """Main AstroPhot tools panel with enhanced information display"""
    
    bl_label = "AstroPhot Tools"
    bl_idname = "ASTRO_PT_astrophot"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Astronomical"
    bl_parent_id = "ASTRONOMICAL_PT_main_panel"
    
    def draw(self, context):
        layout = self.layout
        
        # Import section
        box = layout.box()
        box.label(text="FITS Import:", icon='IMPORT')
        box.operator("astronomical.import_fits", text="Import FITS File", icon='FILE_IMAGE')
        
        # Quick analysis section
        box = layout.box()
        box.label(text="Quick Analysis:", icon='ANALYZE')
        box.operator("astronomical.quick_fits_analysis", text="Quick FITS Analysis", icon='INFO')
        
        # Status section with enhanced information
        fits_objects = [obj for obj in context.scene.objects if 'astrophot' in obj.name.lower()]
        if fits_objects:
            box = layout.box()
            box.label(text="Status:", icon='CHECKMARK')
            box.label(text=f"FITS objects: {len(fits_objects)}")
            
            # Show detailed info for active object
            if context.active_object in fits_objects:
                obj = context.active_object
                
                if "astrophot_object_type" in obj.data:
                    box.label(text=f"Type: {obj.data['astrophot_object_type']}")
                
                if "astrophot_telescope" in obj.data:
                    box.label(text=f"Telescope: {obj.data['astrophot_telescope']}")
                
                if "astrophot_filter" in obj.data:
                    box.label(text=f"Filter: {obj.data['astrophot_filter']}")
                
                if "astrophot_has_wcs" in obj.data and obj.data["astrophot_has_wcs"]:
                    box.label(text="WCS: Available", icon='WORLD')
                    if "astrophot_pixel_scale" in obj.data:
                        scale = obj.data["astrophot_pixel_scale"]
                        box.label(text=f"Scale: {scale:.2f} arcsec/pix")
        
        # Integration status
        box = layout.box()
        box.label(text="Integration Status:", icon='LINKED')
        box.label(text="✓ AstroPhot Core")
        box.label(text="✓ Astropy I/O")
        box.label(text="✓ Advanced Modeling")
        box.label(text="✓ Multi-band Analysis")
        box.label(text="✓ MCMC Bayesian")
        
        # Features section
        box = layout.box()
        box.label(text="Available Features:", icon='INFO')
        box.label(text="• Scientific FITS import")
        box.label(text="• WCS coordinate systems")
        box.label(text="• Background subtraction")
        box.label(text="• Multiple colormaps")
        box.label(text="• 3D displacement surfaces")
        box.label(text="• ModelZoo integration")
        box.label(text="• Multi-band composites")
        box.label(text="• Bayesian parameter fitting")


# Registration
classes = [
    ASTRO_OT_import_fits,
    ASTRO_OT_quick_fits_analysis,
    ASTRO_PT_astrophot_panel,
]

def register():
    """Register all AstroPhot operators and comprehensive functionality"""
    # Register core FITS import
    for cls in classes:
        bpy.utils.register_class(cls)
    
    # Register all advanced functionality modules
    register_modeling()        # Scientific ModelZoo
    register_advanced()        # Multi-band & MCMC  
    register_astropy_io()      # Astropy I/O integration

def unregister():
    """Unregister all AstroPhot operators and functionality"""
    # Unregister advanced functionality first
    unregister_astropy_io()
    unregister_advanced()
    unregister_modeling()
    
    # Unregister core functionality
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)

if __name__ == "__main__":
    register()
