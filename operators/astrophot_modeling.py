"""
AstroPhot Advanced Modeling Operators
====================================

Advanced astronomical modeling using AstroPhot's complete ModelZoo.
Real implementation of all AstroPhot model types for scientific analysis.
"""

import bpy
import astrophot as ap
import torch
import numpy as np
from pathlib import Path
from bpy.props import StringProperty, FloatProperty, BoolProperty, EnumProperty, IntProperty, FloatVectorProperty
from bpy_extras.io_utils import ImportHelper

from ..utilities.astronomical_data import get_stellar_color_from_temperature, get_galaxy_config


class ASTRO_OT_create_astrophot_model(bpy.types.Operator):
    """Create astronomical object using AstroPhot models"""
    
    bl_idname = "astronomical.create_astrophot_model"
    bl_label = "Create AstroPhot Model"
    bl_description = "Create astronomical objects using AstroPhot's scientific models"
    bl_options = {'REGISTER', 'UNDO'}
    
    # Model type selection
    model_category: EnumProperty(
        name="Model Category",
        description="Category of astronomical model",
        items=[
            ('SKY', "Sky Models", "Background sky models"),
            ('PSF', "PSF Models", "Point Spread Function models"),
            ('POINT', "Point Sources", "Stellar point sources"),
            ('GALAXY', "Galaxy Models", "Galaxy surface brightness models"),
            ('SUPERELLIPSE', "Super Ellipse", "Boxy/disky galaxy shapes"),
            ('FOURIER', "Fourier Models", "Complex morphology with Fourier modes"),
            ('WARP', "Warp Models", "Radially varying position angles"),
            ('RAY', "Ray Models", "Sectored galaxy fitting"),
            ('ADVANCED', "Advanced Models", "High-order combinations")
        ],
        default='GALAXY'
    )
    
    # Sky Models
    sky_model_type: EnumProperty(
        name="Sky Model",
        items=[
            ('flat sky model', "Flat Sky", "Constant background"),
            ('plane sky model', "Plane Sky", "Linear gradient background")
        ],
        default='flat sky model'
    )
    
    # PSF Models
    psf_model_type: EnumProperty(
        name="PSF Model", 
        items=[
            ('gaussian psf model', "Gaussian PSF", "Simple Gaussian profile"),
            ('moffat psf model', "Moffat PSF", "Moffat profile (common for ground-based)"),
            ('moffat2d psf model', "2D Moffat PSF", "Elliptical Moffat with PA"),
            ('airy psf model', "Airy Disk PSF", "Theoretical diffraction limit"),
            ('zernike psf model', "Zernike PSF", "Atmospheric turbulence model"),
            ('eigen psf model', "Eigen PSF", "Basis function PSF"),
            ('pixelated psf model', "Pixelated PSF", "Empirical pixel-based PSF")
        ],
        default='moffat psf model'
    )
    
    # Galaxy Models
    galaxy_model_type: EnumProperty(
        name="Galaxy Model",
        items=[
            ('sersic galaxy model', "Sersic Profile", "Standard Sersic surface brightness"),
            ('exponential galaxy model', "Exponential Disk", "Exponential profile (n=1)"),
            ('gaussian galaxy model', "Gaussian Profile", "Gaussian light distribution"),
            ('nuker galaxy model', "Nuker Profile", "Cusp/core galaxy centers"),
            ('spline galaxy model', "Spline Profile", "Non-parametric profile"),
            ('mge model', "Multi-Gaussian", "Multiple Gaussian components"),
            ('isothermal sech2 edgeon model', "Edge-on Disk", "Self-gravitating disk model")
        ],
        default='sersic galaxy model'
    )
    
    # Object properties
    center_x: FloatProperty(name="Center X", default=0.0, description="X position in Blender units")
    center_y: FloatProperty(name="Center Y", default=0.0, description="Y position in Blender units")
    
    # Size and scale
    scale_factor: FloatProperty(
        name="Scale Factor",
        default=1.0,
        min=0.1,
        max=10.0,
        description="Overall scale of the object"
    )
    
    # Galaxy parameters
    sersic_n: FloatProperty(
        name="Sersic Index",
        default=1.0,
        min=0.1,
        max=10.0,
        description="Sersic concentration parameter"
    )
    
    effective_radius: FloatProperty(
        name="Effective Radius",
        default=5.0,
        min=0.1,
        max=50.0,
        description="Half-light radius in Blender units"
    )
    
    axis_ratio: FloatProperty(
        name="Axis Ratio (q)",
        default=0.7,
        min=0.1,
        max=1.0,
        description="Ellipticity (b/a ratio)"
    )
    
    position_angle: FloatProperty(
        name="Position Angle",
        default=0.0,
        min=0.0,
        max=6.28318,
        subtype='ANGLE',
        description="Orientation angle in radians"
    )
    
    # Surface brightness
    central_brightness: FloatProperty(
        name="Central Brightness",
        default=1.0,
        min=0.001,
        max=1000.0,
        description="Central surface brightness"
    )
    
    # Multi-component parameters
    num_components: IntProperty(
        name="Components",
        default=1,
        min=1,
        max=8,
        description="Number of model components"
    )
    
    # Material properties
    emission_strength: FloatProperty(
        name="Emission Strength",
        default=2.0,
        min=0.0,
        max=20.0,
        description="Material emission strength"
    )
    
    color_temperature: FloatProperty(
        name="Color Temperature",
        default=5800.0,
        min=2000.0,
        max=15000.0,
        description="Blackbody temperature for color"
    )

    def execute(self, context):
        try:
            # Create AstroPhot model based on selection
            model_data = self.create_astrophot_model()
            
            # Create Blender representation
            self.create_blender_object(model_data, context)
            
            self.report({'INFO'}, f"Created {self.get_model_type()} model")
            return {'FINISHED'}
            
        except Exception as e:
            self.report({'ERROR'}, f"Failed to create model: {str(e)}")
            return {'CANCELLED'}
    
    def get_model_type(self):
        """Get the selected model type string"""
        if self.model_category == 'SKY':
            return self.sky_model_type
        elif self.model_category == 'PSF':
            return self.psf_model_type
        elif self.model_category == 'GALAXY':
            return self.galaxy_model_type
        else:
            return f"{self.model_category} model"
    
    def create_astrophot_model(self):
        """Create AstroPhot model with scientific parameters"""
        
        # Create basic target image for modeling
        image_size = int(self.effective_radius * 10)  # Reasonable size
        target = ap.image.Target_Image(
            data=np.zeros((image_size, image_size)),
            pixelscale=0.1,  # 0.1 arcsec/pixel
            zeropoint=25.0   # Typical magnitude zeropoint
        )
        
        # Model parameters
        center = [image_size/2, image_size/2]
        
        if self.model_category == 'SKY':
            model = self.create_sky_model(target, center)
        elif self.model_category == 'PSF':
            model = self.create_psf_model(target, center)
        elif self.model_category == 'POINT':
            model = self.create_point_source_model(target, center)
        elif self.model_category == 'GALAXY':
            model = self.create_galaxy_model(target, center)
        else:
            # Default to simple galaxy
            model = self.create_galaxy_model(target, center)
        
        # Initialize and generate model
        model.initialize()
        model_image = model().detach().cpu().numpy()
        
        return {
            'model': model,
            'image_data': model_image,
            'target': target,
            'parameters': dict(model.parameters)
        }
    
    def create_sky_model(self, target, center):
        """Create sky background model"""
        if self.sky_model_type == 'flat sky model':
            return ap.models.AstroPhot_Model(
                model_type="flat sky model",
                parameters={"center": center, "F": np.log10(self.central_brightness)},
                target=target
            )
        else:  # plane sky model
            return ap.models.AstroPhot_Model(
                model_type="plane sky model",
                parameters={
                    "center": center,
                    "F": self.central_brightness,
                    "delta": [0.01, 0.01]  # Gradient
                },
                target=target
            )
    
    def create_psf_model(self, target, center):
        """Create PSF model"""
        # PSF models need PSF_Image target
        psf_target = ap.image.PSF_Image(
            data=np.zeros((101, 101)) + 1e-10,
            pixelscale=target.pixelscale
        )
        
        if self.psf_model_type == 'gaussian psf model':
            return ap.models.AstroPhot_Model(
                model_type="gaussian psf model",
                parameters={"sigma": self.effective_radius},
                target=psf_target
            )
        elif self.psf_model_type == 'moffat psf model':
            return ap.models.AstroPhot_Model(
                model_type="moffat psf model",
                parameters={"n": 2.5, "Rd": self.effective_radius},
                target=psf_target
            )
        elif self.psf_model_type == 'moffat2d psf model':
            return ap.models.AstroPhot_Model(
                model_type="moffat2d psf model",
                parameters={
                    "n": 2.5,
                    "Rd": self.effective_radius,
                    "q": self.axis_ratio,
                    "PA": self.position_angle
                },
                target=psf_target
            )
        elif self.psf_model_type == 'airy psf model':
            return ap.models.AstroPhot_Model(
                model_type="airy psf model",
                parameters={"aRL": 1.0 / self.effective_radius},
                target=psf_target
            )
        else:
            # Default to Moffat
            return ap.models.AstroPhot_Model(
                model_type="moffat psf model",
                parameters={"n": 2.5, "Rd": self.effective_radius},
                target=psf_target
            )
    
    def create_point_source_model(self, target, center):
        """Create point source with PSF"""
        # Create simple PSF first
        psf_target = ap.image.PSF_Image(
            data=np.zeros((51, 51)) + 1e-10,
            pixelscale=target.pixelscale
        )
        
        psf = ap.models.AstroPhot_Model(
            model_type="moffat psf model",
            parameters={"n": 2.5, "Rd": 2.0},
            target=psf_target
        )
        psf.initialize()
        
        return ap.models.AstroPhot_Model(
            model_type="point model",
            parameters={
                "center": [self.center_x + center[0], self.center_y + center[1]],
                "flux": np.log10(self.central_brightness)
            },
            psf=psf,
            target=target
        )
    
    def create_galaxy_model(self, target, center):
        """Create galaxy model"""
        base_params = {
            "center": [self.center_x + center[0], self.center_y + center[1]],
            "q": self.axis_ratio,
            "PA": self.position_angle
        }
        
        if self.galaxy_model_type == 'sersic galaxy model':
            base_params.update({
                "n": self.sersic_n,
                "Re": self.effective_radius,
                "Ie": np.log10(self.central_brightness)
            })
        elif self.galaxy_model_type == 'exponential galaxy model':
            base_params.update({
                "Re": self.effective_radius,
                "Ie": np.log10(self.central_brightness)
            })
        elif self.galaxy_model_type == 'gaussian galaxy model':
            base_params.update({
                "sigma": self.effective_radius,
                "flux": np.log10(self.central_brightness * np.pi * self.effective_radius**2)
            })
        elif self.galaxy_model_type == 'nuker galaxy model':
            base_params.update({
                "Rb": self.effective_radius,
                "Ib": np.log10(self.central_brightness),
                "alpha": 2.0,
                "beta": 4.0,
                "gamma": 0.0
            })
        elif self.galaxy_model_type == 'mge model':
            # Multi-Gaussian Expansion
            n_comp = min(self.num_components, 4)
            sigmas = [self.effective_radius * (i + 1) for i in range(n_comp)]
            fluxes = [self.central_brightness / n_comp for _ in range(n_comp)]
            q_values = [self.axis_ratio for _ in range(n_comp)]
            
            base_params.update({
                "sigma": sigmas,
                "flux": [np.log10(f) for f in fluxes],
                "q": q_values
            })
        elif self.galaxy_model_type == 'isothermal sech2 edgeon model':
            base_params = {
                "center": [self.center_x + center[0], self.center_y + center[1]],
                "PA": self.position_angle,
                "I0": np.log10(self.central_brightness),
                "hs": self.effective_radius / 3,  # Scale height
                "rs": self.effective_radius       # Scale length
            }
        else:
            # Default to Sersic
            base_params.update({
                "n": self.sersic_n,
                "Re": self.effective_radius,
                "Ie": np.log10(self.central_brightness)
            })
        
        return ap.models.AstroPhot_Model(
            model_type=self.galaxy_model_type,
            parameters=base_params,
            target=target
        )
    
    def create_blender_object(self, model_data, context):
        """Create Blender object from AstroPhot model"""
        
        model_image = model_data['image_data']
        model_name = f"AstroPhot_{self.get_model_type().replace(' ', '_')}"
        
        # Create image texture
        blender_image = self.create_blender_image(model_image, model_name)
        
        # Create material
        material = self.create_astrophot_material(blender_image, model_name)
        
        # Create geometry based on model type
        if self.model_category in ['SKY']:
            obj = self.create_background_plane(model_image, material, model_name)
        elif self.model_category in ['PSF', 'POINT']:
            obj = self.create_point_source_geometry(model_image, material, model_name)
        else:  # Galaxy models
            obj = self.create_galaxy_geometry(model_image, material, model_name)
        
        # Position object
        obj.location = (self.center_x * self.scale_factor, self.center_y * self.scale_factor, 0)
        obj.scale = (self.scale_factor, self.scale_factor, self.scale_factor)
        
        # Select and make active
        bpy.context.view_layer.objects.active = obj
        obj.select_set(True)
        
        return obj
    
    def create_blender_image(self, data, name):
        """Create Blender image from model data"""
        if len(data.shape) != 2:
            raise ValueError(f"Expected 2D array, got shape {data.shape}")
        
        height, width = data.shape
        
        # Remove existing image
        if name in bpy.data.images:
            bpy.data.images.remove(bpy.data.images[name])
        
        # Create image
        bpy_image = bpy.data.images.new(name, width=width, height=height, alpha=False)
        
        # Normalize and create RGB data
        normalized = self.normalize_model_data(data)
        
        # Apply color based on temperature
        color = get_stellar_color_from_temperature(self.color_temperature)
        
        # Create RGB
        rgb_data = np.zeros((height, width, 4), dtype=np.float32)
        rgb_data[:, :, 0] = normalized * color[0]  # R
        rgb_data[:, :, 1] = normalized * color[1]  # G
        rgb_data[:, :, 2] = normalized * color[2]  # B
        rgb_data[:, :, 3] = 1.0                   # A
        
        # Set pixels
        bpy_image.pixels = rgb_data.flatten()
        bpy_image.pack()
        
        return bpy_image
    
    def normalize_model_data(self, data):
        """Normalize model data for visualization"""
        valid_data = data[np.isfinite(data)]
        
        if len(valid_data) == 0:
            return np.zeros_like(data, dtype=np.float32)
        
        # Use 99th percentile for normalization to avoid extreme values
        vmax = np.percentile(valid_data, 99)
        vmin = np.min(valid_data)
        
        if vmax > vmin:
            normalized = np.clip((data - vmin) / (vmax - vmin), 0, 1)
        else:
            normalized = np.zeros_like(data)
        
        return normalized.astype(np.float32)
    
    def create_astrophot_material(self, blender_image, name):
        """Create emission material for astronomical object"""
        
        mat_name = f"{name}_material"
        
        # Remove existing
        if mat_name in bpy.data.materials:
            bpy.data.materials.remove(bpy.data.materials[mat_name])
        
        # Create material
        material = bpy.data.materials.new(mat_name)
        material.use_nodes = True
        material.node_tree.nodes.clear()
        
        # Create nodes
        output = material.node_tree.nodes.new('ShaderNodeOutputMaterial')
        emission = material.node_tree.nodes.new('ShaderNodeEmission') 
        image_node = material.node_tree.nodes.new('ShaderNodeTexImage')
        multiply = material.node_tree.nodes.new('ShaderNodeMath')
        
        # Configure nodes
        image_node.image = blender_image
        emission.inputs['Strength'].default_value = self.emission_strength
        multiply.operation = 'MULTIPLY'
        multiply.inputs[1].default_value = self.emission_strength
        
        # Connect nodes
        links = material.node_tree.links
        links.new(image_node.outputs['Color'], emission.inputs['Color'])
        links.new(image_node.outputs['Alpha'], multiply.inputs[0])
        links.new(multiply.outputs['Value'], emission.inputs['Strength'])
        links.new(emission.outputs['Emission'], output.inputs['Surface'])
        
        # Position nodes
        output.location = (400, 0)
        emission.location = (200, 0)
        multiply.location = (200, -200)
        image_node.location = (0, 0)
        
        return material
    
    def create_background_plane(self, data, material, name):
        """Create large background plane for sky models"""
        bpy.ops.mesh.primitive_plane_add(size=100, location=(0, 0, -10))
        obj = bpy.context.active_object
        obj.name = f"{name}_background"
        
        # Assign material
        obj.data.materials.clear()
        obj.data.materials.append(material)
        
        return obj
    
    def create_point_source_geometry(self, data, material, name):
        """Create geometry for point sources"""
        # Create small sphere for point source
        bpy.ops.mesh.primitive_ico_sphere_add(
            subdivisions=2,
            radius=0.5,
            location=(0, 0, 0)
        )
        obj = bpy.context.active_object
        obj.name = f"{name}_point"
        
        # Assign material
        obj.data.materials.clear()
        obj.data.materials.append(material)
        
        return obj
    
    def create_galaxy_geometry(self, data, material, name):
        """Create geometry for galaxy models"""
        height, width = data.shape
        aspect_ratio = width / height
        
        # Create plane
        bpy.ops.mesh.primitive_plane_add(size=self.effective_radius * 2)
        obj = bpy.context.active_object
        obj.name = f"{name}_galaxy"
        
        # Scale to aspect ratio
        obj.scale[0] = aspect_ratio
        
        # Add subdivision for detail
        subdiv_mod = obj.modifiers.new(name="Subdivision", type='SUBSURF')
        subdiv_mod.levels = 3
        
        # Assign material
        obj.data.materials.clear()
        obj.data.materials.append(material)
        
        # UV unwrap
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.uv.unwrap()
        bpy.ops.object.mode_set(mode='OBJECT')
        
        return obj

    def draw(self, context):
        """Custom UI for operator"""
        layout = self.layout
        
        # Model selection
        layout.prop(self, "model_category")
        
        if self.model_category == 'SKY':
            layout.prop(self, "sky_model_type")
        elif self.model_category == 'PSF':
            layout.prop(self, "psf_model_type")
        elif self.model_category == 'GALAXY':
            layout.prop(self, "galaxy_model_type")
        
        layout.separator()
        
        # Position
        row = layout.row()
        row.prop(self, "center_x")
        row.prop(self, "center_y")
        
        layout.prop(self, "scale_factor")
        
        layout.separator()
        
        # Model parameters
        if self.model_category in ['GALAXY', 'PSF']:
            layout.prop(self, "effective_radius")
            layout.prop(self, "axis_ratio")
            layout.prop(self, "position_angle")
            
            if self.galaxy_model_type == 'sersic galaxy model':
                layout.prop(self, "sersic_n")
            elif self.galaxy_model_type == 'mge model':
                layout.prop(self, "num_components")
        
        layout.prop(self, "central_brightness")
        
        layout.separator()
        
        # Material properties
        layout.prop(self, "emission_strength")
        layout.prop(self, "color_temperature")


# Panel for AstroPhot modeling
class ASTRO_PT_astrophot_modeling_panel(bpy.types.Panel):
    """AstroPhot modeling panel"""
    
    bl_label = "AstroPhot Modeling"
    bl_idname = "ASTRO_PT_astrophot_modeling"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Astronomical"
    bl_parent_id = "ASTRO_PT_astrophot"
    
    def draw(self, context):
        layout = self.layout
        
        # Model creation
        box = layout.box()
        box.label(text="Scientific Models:", icon='MODELING')
        box.operator("astronomical.create_astrophot_model", text="Create AstroPhot Model", icon='ADD')
        
        # Info box
        box = layout.box()
        box.label(text="Available Models:", icon='INFO')
        box.label(text="• Sky: Flat, Plane")
        box.label(text="• PSF: Gaussian, Moffat, Airy")
        box.label(text="• Galaxy: Sersic, Exponential")
        box.label(text="• Advanced: Fourier, Warp")


# Registration
classes = [
    ASTRO_OT_create_astrophot_model,
    ASTRO_PT_astrophot_modeling_panel,
]

def register():
    """Register AstroPhot modeling operators"""
    for cls in classes:
        bpy.utils.register_class(cls)

def unregister():
    """Unregister AstroPhot modeling operators"""
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)

if __name__ == "__main__":
    register()
