"""
Stellar Showcase Operator for Astronomical Extension
===================================================

Creates a showcase of different stellar types using modern Blender 4.4 APIs.
No external dependencies - self-contained implementation.
"""

import bpy
from bpy.props import BoolProperty, FloatProperty
from bpy.types import Operator
from mathutils import Vector

class ASTRONOMICAL_OT_CreateStellarShowcase(Operator):
    """Create showcase of different stellar types"""
    bl_idname = "astronomical.create_stellar_showcase"
    bl_label = "Create Stellar Showcase"
    bl_description = "Create showcase of different stellar types with realistic materials"
    bl_options = {"REGISTER", "UNDO"}

    # Layout settings
    spacing: FloatProperty(
        name="Spacing",
        description="Spacing between stars",
        default=4.0,
        min=1.0,
        max=20.0,
    )

    scale: FloatProperty(
        name="Scale",
        description="Scale of stars",
        default=1.0,
        min=0.1,
        max=5.0,
    )

    # Stellar types to include
    include_o_stars: BoolProperty(
        name="O-Type Stars",
        description="Include O-type stars (blue giants)",
        default=True,
    )

    include_b_stars: BoolProperty(
        name="B-Type Stars", 
        description="Include B-type stars (blue-white)",
        default=True,
    )

    include_a_stars: BoolProperty(
        name="A-Type Stars",
        description="Include A-type stars (white)",
        default=True,
    )

    include_f_stars: BoolProperty(
        name="F-Type Stars",
        description="Include F-type stars (yellow-white)",
        default=True,
    )

    include_g_stars: BoolProperty(
        name="G-Type Stars",
        description="Include G-type stars (like our Sun)",
        default=True,
    )

    include_k_stars: BoolProperty(
        name="K-Type Stars",
        description="Include K-type stars (orange)",
        default=True,
    )

    include_m_stars: BoolProperty(
        name="M-Type Stars",
        description="Include M-type stars (red dwarfs)",
        default=True,
    )

    def execute(self, context):
        try:
            # Create stellar showcase
            stars = self._create_stellar_showcase()

            self.report({"INFO"}, f"Created stellar showcase with {len(stars)} stars!")
            return {"FINISHED"}

        except Exception as e:
            self.report({"ERROR"}, f"Failed to create stellar showcase: {str(e)}")
            return {"CANCELLED"}

    def _create_stellar_showcase(self):
        """Create showcase of different stellar types"""
        stellar_types = []
        x_offset = 0

        # Define stellar types with realistic properties
        stellar_data = {
            'O': {'temp': 30000, 'color': (0.6, 0.7, 1.0), 'size': 2.0, 'brightness': 10.0},
            'B': {'temp': 20000, 'color': (0.7, 0.8, 1.0), 'size': 1.8, 'brightness': 8.0},
            'A': {'temp': 8500, 'color': (0.9, 0.9, 1.0), 'size': 1.5, 'brightness': 6.0},
            'F': {'temp': 6500, 'color': (1.0, 1.0, 0.9), 'size': 1.3, 'brightness': 4.0},
            'G': {'temp': 5500, 'color': (1.0, 1.0, 0.8), 'size': 1.0, 'brightness': 3.0},  # Sun
            'K': {'temp': 4000, 'color': (1.0, 0.8, 0.6), 'size': 0.8, 'brightness': 2.0},
            'M': {'temp': 3000, 'color': (1.0, 0.5, 0.3), 'size': 0.5, 'brightness': 1.0},
        }

        # Build list based on user selection
        if self.include_o_stars:
            stellar_types.append(('O', x_offset))
            x_offset += self.spacing

        if self.include_b_stars:
            stellar_types.append(('B', x_offset))
            x_offset += self.spacing

        if self.include_a_stars:
            stellar_types.append(('A', x_offset))
            x_offset += self.spacing

        if self.include_f_stars:
            stellar_types.append(('F', x_offset))
            x_offset += self.spacing

        if self.include_g_stars:
            stellar_types.append(('G', x_offset))
            x_offset += self.spacing

        if self.include_k_stars:
            stellar_types.append(('K', x_offset))
            x_offset += self.spacing

        if self.include_m_stars:
            stellar_types.append(('M', x_offset))
            x_offset += self.spacing

        if not stellar_types:
            self.report({"WARNING"}, "No stellar types selected!")
            return []

        # Create stars
        stars = []
        for star_class, x_pos in stellar_types:
            star_obj = self._create_single_star(star_class, stellar_data[star_class], x_pos)
            stars.append(star_obj)

        # Create text labels
        self._create_star_labels(stellar_types, stellar_data)

        return stars

    def _create_single_star(self, star_class, properties, x_position):
        """Create a single star with realistic properties"""
        # Create sphere for star
        bpy.ops.mesh.primitive_uv_sphere_add(
            radius=properties['size'] * self.scale,
            location=(x_position, 0, 0)
        )
        
        star_obj = bpy.context.active_object
        star_obj.name = f"Star_{star_class}_Type"
        
        # Create material
        mat = self._create_star_material(star_class, properties)
        star_obj.data.materials.append(mat)
        
        # Add emission for glow
        self._add_star_glow(star_obj, properties)
        
        return star_obj

    def _create_star_material(self, star_class, properties):
        """Create realistic star material"""
        mat_name = f"Star_Material_{star_class}"
        mat = bpy.data.materials.new(mat_name)
        mat.use_nodes = True
        
        # Clear default nodes
        mat.node_tree.nodes.clear()
        
        # Create emission shader
        emission = mat.node_tree.nodes.new('ShaderNodeEmission')
        emission.location = (0, 0)
        emission.inputs['Color'].default_value = (*properties['color'], 1.0)
        emission.inputs['Strength'].default_value = properties['brightness']
        
        # Create output
        output = mat.node_tree.nodes.new('ShaderNodeOutputMaterial')
        output.location = (200, 0)
        
        # Link nodes
        mat.node_tree.links.new(emission.outputs['Emission'], output.inputs['Surface'])
        
        # Optional: Add blackbody node for realistic color
        if bpy.app.version >= (3, 0, 0):
            blackbody = mat.node_tree.nodes.new('ShaderNodeBlackbody')
            blackbody.location = (-200, 0)
            blackbody.inputs['Temperature'].default_value = properties['temp']
            
            # Mix with emission color
            mix = mat.node_tree.nodes.new('ShaderNodeMix')
            mix.location = (-100, 0)
            mix.data_type = 'RGBA'
            mix.inputs['Fac'].default_value = 0.8  # 80% blackbody, 20% artistic
            
            mat.node_tree.links.new(blackbody.outputs['Color'], mix.inputs['Color1'])
            mat.node_tree.links.new(mix.outputs['Result'], emission.inputs['Color'])
        
        return mat

    def _add_star_glow(self, star_obj, properties):
        """Add glow effect around star"""
        # Create larger transparent sphere for glow
        bpy.ops.mesh.primitive_uv_sphere_add(
            radius=properties['size'] * self.scale * 2.0,
            location=star_obj.location
        )
        
        glow_obj = bpy.context.active_object
        glow_obj.name = f"{star_obj.name}_Glow"
        
        # Create glow material
        glow_mat = bpy.data.materials.new(f"{star_obj.name}_Glow_Material")
        glow_mat.use_nodes = True
        glow_mat.blend_method = 'ADD'  # Additive blending
        
        nodes = glow_mat.node_tree.nodes
        nodes.clear()
        
        # Emission for glow
        emission = nodes.new('ShaderNodeEmission')
        emission.inputs['Color'].default_value = (*properties['color'], 1.0)
        emission.inputs['Strength'].default_value = properties['brightness'] * 0.2
        
        # Fresnel for edge glow
        fresnel = nodes.new('ShaderNodeFresnel')
        fresnel.inputs['IOR'].default_value = 1.1
        
        # Mix with fresnel
        mix = nodes.new('ShaderNodeMix')
        mix.data_type = 'RGBA'
        
        # Transparent shader
        transparent = nodes.new('ShaderNodeBsdfTransparent')
        
        # Mix shader
        mix_shader = nodes.new('ShaderNodeMixShader')
        
        # Output
        output = nodes.new('ShaderNodeOutputMaterial')
        
        # Links
        glow_mat.node_tree.links.new(fresnel.outputs['Fac'], mix_shader.inputs['Fac'])
        glow_mat.node_tree.links.new(transparent.outputs['BSDF'], mix_shader.inputs[1])
        glow_mat.node_tree.links.new(emission.outputs['Emission'], mix_shader.inputs[2])
        glow_mat.node_tree.links.new(mix_shader.outputs['Shader'], output.inputs['Surface'])
        
        glow_obj.data.materials.append(glow_mat)
        
        # Parent glow to star
        glow_obj.parent = star_obj

    def _create_star_labels(self, stellar_types, stellar_data):
        """Create text labels for each star type"""
        for star_class, x_pos in stellar_types:
            # Create text object
            bpy.ops.object.text_add(location=(x_pos, -2, 0))
            text_obj = bpy.context.active_object
            text_obj.name = f"Label_{star_class}_Type"
            
            # Set text content
            text_obj.data.body = f"{star_class}-Type\n{stellar_data[star_class]['temp']}K"
            text_obj.data.align_x = 'CENTER'
            text_obj.data.align_y = 'CENTER'
            
            # Scale down text
            text_obj.scale = (0.5, 0.5, 0.5)
            
            # Create text material
            text_mat = bpy.data.materials.new(f"Label_Material_{star_class}")
            text_mat.use_nodes = True
            
            # Simple emission for visibility
            emission = text_mat.node_tree.nodes["Emission"]
            emission.inputs['Color'].default_value = (1.0, 1.0, 1.0, 1.0)
            emission.inputs['Strength'].default_value = 2.0
            
            text_obj.data.materials.append(text_mat)


def register():
    bpy.utils.register_class(ASTRONOMICAL_OT_CreateStellarShowcase)


def unregister():
    bpy.utils.unregister_class(ASTRONOMICAL_OT_CreateStellarShowcase)
