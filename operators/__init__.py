"""
Astronomical Extension Operators
===============================

Modern Blender 4.4 operators following clean refactoring principles.
Includes AstroPhot integration for FITS data processing.
"""

import bpy
from bpy.types import Operator
from bpy.props import IntProperty, FloatProperty, EnumProperty, BoolProperty

# Import AstroPhot operators
from .astrophot import register as register_astrophot
from .astrophot import unregister as unregister_astrophot

# ==========================================
# Core Astronomical Operators  
# ==========================================

class ASTRONOMICAL_OT_create_star_field(Operator):
    """Create a field of stars with realistic distribution"""
    bl_idname = "astronomical.create_star_field"
    bl_label = "Create Star Field"
    bl_description = "Generate a field of stars using modern geometry nodes"
    bl_options = {'REGISTER', 'UNDO'}

    count: IntProperty(
        name="Star Count",
        description="Number of stars to generate",
        default=1000,
        min=10,
        max=100000
    )

    distribution: EnumProperty(
        name="Distribution",
        description="Star distribution pattern",
        items=[
            ('RANDOM', "Random", "Random distribution"),
            ('CLUSTERED', "Clustered", "Clustered around center"),
            ('SPIRAL', "Spiral", "Spiral galaxy pattern"),
            ('DISK', "Disk", "Flat disk distribution"),
        ],
        default='RANDOM'
    )

    scale: FloatProperty(
        name="Scale",
        description="Scale factor for star sizes",
        default=1.0,
        min=0.1,
        max=10.0
    )

    radius: FloatProperty(
        name="Field Radius",
        description="Radius of star field",
        default=50.0,
        min=1.0,
        max=1000.0
    )

    def execute(self, context):
        """Execute star field creation"""
        try:
            # Create star field using modern Blender 4.4 approach
            self._create_star_field_geometry_nodes(context)
            
            self.report({'INFO'}, f"Created star field with {self.count} stars")
            return {'FINISHED'}
            
        except Exception as e:
            self.report({'ERROR'}, f"Failed to create star field: {e}")
            return {'CANCELLED'}

    def _create_star_field_geometry_nodes(self, context):
        """Create star field using geometry nodes (modern approach)"""
        import bmesh
        from mathutils import Vector
        import random
        
        # Create base mesh
        mesh = bpy.data.meshes.new("StarField")
        obj = bpy.data.objects.new("StarField", mesh)
        context.collection.objects.link(obj)
        
        # Create geometry with bmesh
        bm = bmesh.new()
        
        for i in range(self.count):
            if self.distribution == 'RANDOM':
                # Random spherical distribution
                pos = Vector((
                    random.uniform(-self.radius, self.radius),
                    random.uniform(-self.radius, self.radius), 
                    random.uniform(-self.radius, self.radius)
                ))
            elif self.distribution == 'DISK':
                # Disk distribution (galactic)
                import math
                angle = random.uniform(0, 6.28318)  # 2*pi
                r = random.uniform(0, self.radius)
                pos = Vector((
                    r * random.uniform(-1, 1) * 0.1,  # Thin disk
                    r * math.cos(angle),
                    r * math.sin(angle)
                ))
            else:
                # Default to random
                pos = Vector((
                    random.uniform(-self.radius, self.radius),
                    random.uniform(-self.radius, self.radius),
                    random.uniform(-self.radius, self.radius)
                ))
            
            # Create vertex for each star
            bm.verts.new(pos)
        
        # Update mesh
        bm.to_mesh(mesh)
        bm.free()
        
        # Add geometry nodes modifier for instancing
        self._add_star_geometry_nodes(obj)

    def _add_star_geometry_nodes(self, obj):
        """Add geometry nodes for star instancing"""
        # Add geometry nodes modifier
        modifier = obj.modifiers.new("StarInstancing", 'NODES')
        
        # Create node group
        node_group = bpy.data.node_groups.new("StarField", 'GeometryNodeTree')
        modifier.node_group = node_group
        
        # Create nodes using modern Blender 4.4 Interface API
        nodes = node_group.nodes
        links = node_group.links
        
        # Input node
        input_node = nodes.new('NodeGroupInput')
        input_node.location = (-400, 0)
        
        # Output node  
        output_node = nodes.new('NodeGroupOutput')
        output_node.location = (400, 0)
        
        # Create sockets using modern Interface API
        node_group.interface.new_socket(
            name='Geometry', 
            in_out='INPUT', 
            socket_type='NodeSocketGeometry'
        )
        node_group.interface.new_socket(
            name='Geometry',
            in_out='OUTPUT', 
            socket_type='NodeSocketGeometry'
        )
        
        # Instance on Points node
        instance_node = nodes.new('GeometryNodeInstanceOnPoints')
        instance_node.location = (0, 0)
        
        # Create star mesh to instance
        star_mesh = self._create_star_mesh()
        instance_node.inputs['Instance'].default_value = star_mesh
        
        # Link nodes
        links.new(input_node.outputs['Geometry'], instance_node.inputs['Points'])
        links.new(instance_node.outputs['Instances'], output_node.inputs['Geometry'])

    def _create_star_mesh(self):
        """Create basic star mesh for instancing"""
        # Create simple star mesh
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=1, radius=0.1)
        star_obj = bpy.context.active_object
        star_obj.name = "StarInstance"
        
        # Add emission material
        mat = bpy.data.materials.new("StarMaterial")
        mat.use_nodes = True
        
        # Clear default nodes
        mat.node_tree.nodes.clear()
        
        # Create emission shader
        emission = mat.node_tree.nodes.new('ShaderNodeEmission')
        emission.inputs['Color'].default_value = (1.0, 1.0, 0.8, 1.0)  # Warm white
        emission.inputs['Strength'].default_value = 5.0
        
        # Create output
        output = mat.node_tree.nodes.new('ShaderNodeOutputMaterial')
        
        # Link
        mat.node_tree.links.new(emission.outputs['Emission'], output.inputs['Surface'])
        
        # Assign material
        star_obj.data.materials.append(mat)
        
        return star_obj

class ASTRONOMICAL_OT_create_galaxy(Operator):
    """Create a galaxy using modern geometry nodes"""
    bl_idname = "astronomical.create_galaxy"
    bl_label = "Create Galaxy"
    bl_description = "Generate a galaxy with spiral arms or elliptical shape"
    bl_options = {'REGISTER', 'UNDO'}

    galaxy_type: EnumProperty(
        name="Galaxy Type",
        description="Type of galaxy to create",
        items=[
            ('SPIRAL', "Spiral", "Spiral galaxy with arms"),
            ('ELLIPTICAL', "Elliptical", "Elliptical galaxy"),
            ('IRREGULAR', "Irregular", "Irregular galaxy"),
        ],
        default='SPIRAL'
    )

    arms: IntProperty(
        name="Spiral Arms",
        description="Number of spiral arms",
        default=4,
        min=2,
        max=8
    )

    bulge_size: FloatProperty(
        name="Bulge Size", 
        description="Size of galactic bulge",
        default=1.0,
        min=0.1,
        max=5.0
    )

    def execute(self, context):
        """Execute galaxy creation"""
        try:
            self._create_galaxy(context)
            self.report({'INFO'}, f"Created {self.galaxy_type.lower()} galaxy")
            return {'FINISHED'}
            
        except Exception as e:
            self.report({'ERROR'}, f"Failed to create galaxy: {e}")
            return {'CANCELLED'}

    def _create_galaxy(self, context):
        """Create galaxy geometry"""
        # Implementation would use geometry nodes for galaxy generation
        # For now, create placeholder
        bpy.ops.mesh.primitive_uv_sphere_add(radius=10)
        galaxy_obj = context.active_object
        galaxy_obj.name = f"Galaxy_{self.galaxy_type}"
        
        # Add galaxy material
        self._add_galaxy_material(galaxy_obj)

    def _add_galaxy_material(self, obj):
        """Add galaxy material with emission"""
        mat = bpy.data.materials.new("GalaxyMaterial")
        mat.use_nodes = True
        
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links
        
        # Clear default
        nodes.clear()
        
        # Create galaxy shader
        emission = nodes.new('ShaderNodeEmission')
        emission.inputs['Color'].default_value = (0.8, 0.6, 1.0, 1.0)  # Purple-blue
        emission.inputs['Strength'].default_value = 2.0
        
        output = nodes.new('ShaderNodeOutputMaterial')
        links.new(emission.outputs['Emission'], output.inputs['Surface'])
        
        obj.data.materials.append(mat)

class ASTRONOMICAL_OT_setup_scene(Operator):
    """Setup scene for astronomical rendering"""
    bl_idname = "astronomical.setup_scene"
    bl_label = "Setup Astronomical Scene"
    bl_description = "Configure scene for astronomical visualization"
    bl_options = {'REGISTER', 'UNDO'}

    preset: EnumProperty(
        name="Preset",
        description="Scene preset",
        items=[
            ('SCIENTIFIC', "Scientific", "Black background, minimal lighting"),
            ('CINEMATIC', "Cinematic", "Dramatic lighting and effects"),
            ('PUBLICATION', "Publication", "White background for papers"),
        ],
        default='SCIENTIFIC'
    )

    def execute(self, context):
        """Execute scene setup"""
        try:
            from ..core.scene import setup_scene
            
            # Map presets
            preset_map = {
                'SCIENTIFIC': 'scientific',
                'CINEMATIC': 'cinematic', 
                'PUBLICATION': 'publication'
            }
            
            setup_scene(preset=preset_map[self.preset])
            
            # Setup render engine for Blender 4.4
            context.scene.render.engine = 'BLENDER_EEVEE_NEXT'
            
            self.report({'INFO'}, f"Scene configured for {self.preset.lower()} use")
            return {'FINISHED'}
            
        except Exception as e:
            self.report({'ERROR'}, f"Failed to setup scene: {e}")
            return {'CANCELLED'}

# ==========================================
# UI Panel
# ==========================================

class ASTRONOMICAL_PT_main_panel(bpy.types.Panel):
    """Main panel for astronomical tools"""
    bl_label = "Astronomical Tools"
    bl_idname = "ASTRONOMICAL_PT_main_panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Astronomical"

    def draw(self, context):
        layout = self.layout
        
        # Scene setup
        box = layout.box()
        box.label(text="Scene Setup", icon='SCENE_DATA')
        box.operator("astronomical.setup_scene")
        
        layout.separator()
        
        # Object creation
        box = layout.box()
        box.label(text="Create Objects", icon='MESH_DATA')
        box.operator("astronomical.create_star_field", icon='ADD')
        box.operator("astronomical.create_galaxy", icon='MESH_UVSPHERE')
        
        layout.separator()
        
        # Data import section (AstroPhot integration)
        box = layout.box()
        box.label(text="Scientific Data", icon='IMPORT')
        box.operator("astronomical.import_fits", text="Import FITS", icon='FILE_IMAGE')
        
        layout.separator()
        
        # Render info
        box = layout.box()
        box.label(text="Render Engine", icon='RENDER_STILL')
        scene = context.scene
        box.label(text=f"Engine: {scene.render.engine}")
        
        if scene.render.engine == 'BLENDER_EEVEE_NEXT':
            box.label(text="✓ Modern EEVEE Next", icon='CHECKMARK')
        else:
            box.label(text="⚠ Consider EEVEE Next", icon='ERROR')

# ==========================================
# Registration
# ==========================================

classes = [
    ASTRONOMICAL_OT_create_star_field,
    ASTRONOMICAL_OT_create_galaxy,
    ASTRONOMICAL_OT_setup_scene,
    ASTRONOMICAL_PT_main_panel,
]

def register():
    """Register all operators"""
    print("Registering Astronomical Operators...")
    
    # Register core operators
    for cls in classes:
        bpy.utils.register_class(cls)
    
    # Register AstroPhot integration
    register_astrophot()
    
    print("✅ Astronomical Operators registered (including AstroPhot)")

def unregister():
    """Unregister all operators"""
    print("Unregistering Astronomical Operators...")
    
    # Unregister AstroPhot integration
    unregister_astrophot()
    
    # Unregister core operators
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)
    
    print("✅ Astronomical Operators unregistered")

if __name__ == "__main__":
    register()
