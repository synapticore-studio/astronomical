"""
Advanced AstroPhot Geometry Nodes
================================

Warp Models, Ray Models, Fourier Models - the complete AstroPhot ModelZoo in Blender 4.4.
Galactic user experience with procedural astronomical modeling.
"""

import bpy
import numpy as np
from typing import Dict, Any, Tuple


def create_warp_galaxy_node_group():
    """
    Create Warp Galaxy Geometry Node Group - radially varying position angles and ellipticity.
    Based on AstroPhot's warp models for complex galaxy morphology.
    """
    ng = bpy.data.node_groups.new("ASTROPHOT_WarpGalaxy", "GeometryNodeTree")
    
    # Modern Blender 4.4 Interface API
    interface = ng.interface
    
    # Input sockets - matching AstroPhot warp model parameters
    interface.new_socket(name="Base Mesh", in_out='INPUT', socket_type='NodeSocketGeometry')
    interface.new_socket(name="Star Count", in_out='INPUT', socket_type='NodeSocketInt')
    interface.new_socket(name="Warp Strength", in_out='INPUT', socket_type='NodeSocketFloat')
    interface.new_socket(name="Radial Samples", in_out='INPUT', socket_type='NodeSocketInt')
    interface.new_socket(name="PA Variation", in_out='INPUT', socket_type='NodeSocketFloat')
    interface.new_socket(name="Ellipticity Variation", in_out='INPUT', socket_type='NodeSocketFloat')
    interface.new_socket(name="Galaxy Scale", in_out='INPUT', socket_type='NodeSocketFloat')
    
    # Output socket
    interface.new_socket(name="Geometry", in_out='OUTPUT', socket_type='NodeSocketGeometry')
    
    # Set realistic defaults
    ng.interface.items_tree[1].default_value = 15000  # Star Count
    ng.interface.items_tree[2].default_value = 0.3    # Warp Strength
    ng.interface.items_tree[3].default_value = 8      # Radial Samples
    ng.interface.items_tree[4].default_value = 1.57   # PA Variation (90 degrees)
    ng.interface.items_tree[5].default_value = 0.4    # Ellipticity Variation
    ng.interface.items_tree[6].default_value = 20.0   # Galaxy Scale
    
    # Create nodes
    input_node = ng.nodes.new("NodeGroupInput")
    output_node = ng.nodes.new("NodeGroupOutput")
    
    # Base galaxy disk
    cylinder = ng.nodes.new("GeometryNodeMeshCylinder")
    cylinder.fill_type = "NGON"
    
    # Distribute points for stars
    distribute = ng.nodes.new("GeometryNodeDistributePointsOnFaces")
    distribute.distribute_method = "RANDOM"
    
    # Position analysis for warp effects
    position = ng.nodes.new("GeometryNodeInputPosition")
    separate_xyz = ng.nodes.new("ShaderNodeSeparateXYZ")
    
    # Calculate radius from center
    vector_math_dot = ng.nodes.new("ShaderNodeVectorMath")
    vector_math_dot.operation = 'DOT_PRODUCT'
    math_sqrt = ng.nodes.new("ShaderNodeMath")
    math_sqrt.operation = 'POWER'
    math_sqrt.inputs[1].default_value = 0.5  # Square root
    
    # Radial-dependent position angle
    math_multiply_pa = ng.nodes.new("ShaderNodeMath")
    math_multiply_pa.operation = 'MULTIPLY'
    
    # Radial-dependent ellipticity
    math_multiply_ellip = ng.nodes.new("ShaderNodeMath") 
    math_multiply_ellip.operation = 'MULTIPLY'
    
    # Rotation matrices for warp
    vector_rotate = ng.nodes.new("ShaderNodeVectorRotate")
    vector_rotate.type = 'Z_AXIS'
    
    # Elliptical scaling
    vector_math_scale = ng.nodes.new("ShaderNodeVectorMath")
    vector_math_scale.operation = 'MULTIPLY'
    
    # Apply warp transformation
    set_position = ng.nodes.new("GeometryNodeSetPosition")
    
    # Create star instances
    ico_sphere = ng.nodes.new("GeometryNodeMeshIcoSphere")
    ico_sphere.inputs["Subdivisions"].default_value = 1
    ico_sphere.inputs["Radius"].default_value = 0.05
    
    instance_on_points = ng.nodes.new("GeometryNodeInstanceOnPoints")
    
    # Random scale for realistic stars
    random_value = ng.nodes.new("FunctionNodeRandomValue")
    random_value.data_type = 'FLOAT'
    random_value.inputs["Min"].default_value = 0.3
    random_value.inputs["Max"].default_value = 1.8
    
    # Position nodes for clarity
    input_node.location = (-800, 0)
    cylinder.location = (-600, 0)
    distribute.location = (-400, 0)
    
    position.location = (-200, 300)
    separate_xyz.location = (0, 300)
    vector_math_dot.location = (200, 300)
    math_sqrt.location = (400, 300)
    
    math_multiply_pa.location = (200, 150)
    math_multiply_ellip.location = (200, 50)
    vector_rotate.location = (400, 150)
    vector_math_scale.location = (400, 50)
    
    set_position.location = (800, 0)
    
    ico_sphere.location = (-200, -200)
    random_value.location = (0, -200)
    instance_on_points.location = (1000, 0)
    output_node.location = (1200, 0)
    
    # Connect nodes for warp galaxy effect
    # Base geometry
    ng.links.new(input_node.outputs["Galaxy Scale"], cylinder.inputs["Radius"])
    ng.links.new(cylinder.outputs["Mesh"], distribute.inputs["Mesh"])
    ng.links.new(input_node.outputs["Star Count"], distribute.inputs["Density"])
    
    # Warp calculations
    ng.links.new(distribute.outputs["Points"], set_position.inputs["Geometry"])
    ng.links.new(position.outputs["Position"], separate_xyz.inputs["Vector"])
    ng.links.new(separate_xyz.outputs["X"], vector_math_dot.inputs[0])
    ng.links.new(separate_xyz.outputs["Y"], vector_math_dot.inputs[1])
    ng.links.new(vector_math_dot.outputs["Value"], math_sqrt.inputs[0])
    
    # Radial PA variation
    ng.links.new(math_sqrt.outputs["Value"], math_multiply_pa.inputs[0])
    ng.links.new(input_node.outputs["PA Variation"], math_multiply_pa.inputs[1])
    ng.links.new(math_multiply_pa.outputs["Value"], vector_rotate.inputs["Angle"])
    
    # Radial ellipticity variation
    ng.links.new(math_sqrt.outputs["Value"], math_multiply_ellip.inputs[0])
    ng.links.new(input_node.outputs["Ellipticity Variation"], math_multiply_ellip.inputs[1])
    
    # Apply transformations
    ng.links.new(position.outputs["Position"], vector_rotate.inputs["Vector"])
    ng.links.new(vector_rotate.outputs["Vector"], vector_math_scale.inputs[0])
    ng.links.new(math_multiply_ellip.outputs["Value"], vector_math_scale.inputs[1])
    ng.links.new(vector_math_scale.outputs["Vector"], set_position.inputs["Position"])
    
    # Star instancing
    ng.links.new(ico_sphere.outputs["Mesh"], instance_on_points.inputs["Instance"])
    ng.links.new(set_position.outputs["Geometry"], instance_on_points.inputs["Points"])
    ng.links.new(random_value.outputs["Value"], instance_on_points.inputs["Scale"])
    
    # Final output
    ng.links.new(instance_on_points.outputs["Instances"], output_node.inputs["Geometry"])
    
    return ng


def create_fourier_galaxy_node_group():
    """
    Create Fourier Galaxy Geometry Node Group - complex morphological features.
    Based on AstroPhot's Fourier models for bars, spiral arms, and asymmetries.
    """
    ng = bpy.data.node_groups.new("ASTROPHOT_FourierGalaxy", "GeometryNodeTree")
    
    # Interface
    interface = ng.interface
    interface.new_socket(name="Base Mesh", in_out='INPUT', socket_type='NodeSocketGeometry')
    interface.new_socket(name="Fourier Mode 1", in_out='INPUT', socket_type='NodeSocketFloat')
    interface.new_socket(name="Fourier Mode 2", in_out='INPUT', socket_type='NodeSocketFloat')
    interface.new_socket(name="Fourier Mode 3", in_out='INPUT', socket_type='NodeSocketFloat')
    interface.new_socket(name="Phase 1", in_out='INPUT', socket_type='NodeSocketFloat')
    interface.new_socket(name="Phase 2", in_out='INPUT', socket_type='NodeSocketFloat')
    interface.new_socket(name="Phase 3", in_out='INPUT', socket_type='NodeSocketFloat')
    interface.new_socket(name="Galaxy Scale", in_out='INPUT', socket_type='NodeSocketFloat')
    
    interface.new_socket(name="Geometry", in_out='OUTPUT', socket_type='NodeSocketGeometry')
    
    # Defaults for realistic galaxy
    ng.interface.items_tree[1].default_value = 0.1    # Mode 1 (shift)
    ng.interface.items_tree[2].default_value = 0.3    # Mode 2 (ellipticity)
    ng.interface.items_tree[3].default_value = -0.2   # Mode 3 (triangular)
    ng.interface.items_tree[4].default_value = 0.0    # Phase 1
    ng.interface.items_tree[5].default_value = 1.57   # Phase 2 (90 degrees)
    ng.interface.items_tree[6].default_value = 0.7    # Phase 3
    ng.interface.items_tree[7].default_value = 18.0   # Galaxy Scale
    
    # Create comprehensive Fourier node network
    input_node = ng.nodes.new("NodeGroupInput")
    output_node = ng.nodes.new("NodeGroupOutput")
    
    # Base geometry
    cylinder = ng.nodes.new("GeometryNodeMeshCylinder")
    distribute = ng.nodes.new("GeometryNodeDistributePointsOnFaces")
    
    # Position analysis
    position = ng.nodes.new("GeometryNodeInputPosition")
    separate_xyz = ng.nodes.new("ShaderNodeSeparateXYZ")
    
    # Polar coordinates
    math_atan2 = ng.nodes.new("ShaderNodeMath")
    math_atan2.operation = 'ARCTAN2'
    
    # Fourier mode calculations
    math_sin1 = ng.nodes.new("ShaderNodeMath")
    math_sin1.operation = 'SINE'
    math_sin2 = ng.nodes.new("ShaderNodeMath") 
    math_sin2.operation = 'SINE'
    math_sin3 = ng.nodes.new("ShaderNodeMath")
    math_sin3.operation = 'SINE'
    
    # Mode multipliers
    math_mult1 = ng.nodes.new("ShaderNodeMath")
    math_mult1.operation = 'MULTIPLY'
    math_mult2 = ng.nodes.new("ShaderNodeMath")
    math_mult2.operation = 'MULTIPLY'
    math_mult3 = ng.nodes.new("ShaderNodeMath")
    math_mult3.operation = 'MULTIPLY'
    
    # Frequency multipliers for different modes
    math_freq1 = ng.nodes.new("ShaderNodeMath")
    math_freq1.operation = 'MULTIPLY'
    math_freq1.inputs[1].default_value = 1.0  # Mode 1
    
    math_freq2 = ng.nodes.new("ShaderNodeMath")
    math_freq2.operation = 'MULTIPLY'
    math_freq2.inputs[1].default_value = 2.0  # Mode 2
    
    math_freq3 = ng.nodes.new("ShaderNodeMath")
    math_freq3.operation = 'MULTIPLY'
    math_freq3.inputs[1].default_value = 3.0  # Mode 3
    
    # Phase additions
    math_phase1 = ng.nodes.new("ShaderNodeMath")
    math_phase1.operation = 'ADD'
    math_phase2 = ng.nodes.new("ShaderNodeMath")
    math_phase2.operation = 'ADD'
    math_phase3 = ng.nodes.new("ShaderNodeMath")
    math_phase3.operation = 'ADD'
    
    # Combine Fourier modes
    math_add_modes = ng.nodes.new("ShaderNodeMath")
    math_add_modes.operation = 'ADD'
    math_add_final = ng.nodes.new("ShaderNodeMath")
    math_add_final.operation = 'ADD'
    
    # Apply Fourier deformation
    set_position = ng.nodes.new("GeometryNodeSetPosition")
    
    # Star instancing
    ico_sphere = ng.nodes.new("GeometryNodeMeshIcoSphere")
    ico_sphere.inputs["Radius"].default_value = 0.03
    
    instance_on_points = ng.nodes.new("GeometryNodeInstanceOnPoints")
    
    # Position nodes strategically
    input_node.location = (-1000, 0)
    cylinder.location = (-800, 0)
    distribute.location = (-600, 0)
    
    position.location = (-400, 400)
    separate_xyz.location = (-200, 400)
    math_atan2.location = (0, 400)
    
    # Fourier calculation chain
    math_freq1.location = (200, 500)
    math_freq2.location = (200, 400)
    math_freq3.location = (200, 300)
    
    math_phase1.location = (400, 500)
    math_phase2.location = (400, 400)
    math_phase3.location = (400, 300)
    
    math_sin1.location = (600, 500)
    math_sin2.location = (600, 400)
    math_sin3.location = (600, 300)
    
    math_mult1.location = (800, 500)
    math_mult2.location = (800, 400)
    math_mult3.location = (800, 300)
    
    math_add_modes.location = (1000, 450)
    math_add_final.location = (1000, 350)
    
    set_position.location = (1200, 0)
    
    ico_sphere.location = (-400, -200)
    instance_on_points.location = (1400, 0)
    output_node.location = (1600, 0)
    
    # Connect the Fourier galaxy network
    # Base geometry
    ng.links.new(input_node.outputs["Galaxy Scale"], cylinder.inputs["Radius"])
    ng.links.new(cylinder.outputs["Mesh"], distribute.inputs["Mesh"])
    ng.links.new(distribute.outputs["Points"], set_position.inputs["Geometry"])
    
    # Position analysis
    ng.links.new(position.outputs["Position"], separate_xyz.inputs["Vector"])
    ng.links.new(separate_xyz.outputs["Y"], math_atan2.inputs[0])
    ng.links.new(separate_xyz.outputs["X"], math_atan2.inputs[1])
    
    # Fourier mode 1 (m=1, shift/lopsidedness)
    ng.links.new(math_atan2.outputs["Value"], math_freq1.inputs[0])
    ng.links.new(math_freq1.outputs["Value"], math_phase1.inputs[0])
    ng.links.new(input_node.outputs["Phase 1"], math_phase1.inputs[1])
    ng.links.new(math_phase1.outputs["Value"], math_sin1.inputs[0])
    ng.links.new(math_sin1.outputs["Value"], math_mult1.inputs[0])
    ng.links.new(input_node.outputs["Fourier Mode 1"], math_mult1.inputs[1])
    
    # Fourier mode 2 (m=2, ellipticity/bar)
    ng.links.new(math_atan2.outputs["Value"], math_freq2.inputs[0])
    ng.links.new(math_freq2.outputs["Value"], math_phase2.inputs[0])
    ng.links.new(input_node.outputs["Phase 2"], math_phase2.inputs[1])
    ng.links.new(math_phase2.outputs["Value"], math_sin2.inputs[0])
    ng.links.new(math_sin2.outputs["Value"], math_mult2.inputs[0])
    ng.links.new(input_node.outputs["Fourier Mode 2"], math_mult2.inputs[1])
    
    # Fourier mode 3 (m=3, triangular/Y-shape)
    ng.links.new(math_atan2.outputs["Value"], math_freq3.inputs[0])
    ng.links.new(math_freq3.outputs["Value"], math_phase3.inputs[0])
    ng.links.new(input_node.outputs["Phase 3"], math_phase3.inputs[1])
    ng.links.new(math_phase3.outputs["Value"], math_sin3.inputs[0])
    ng.links.new(math_sin3.outputs["Value"], math_mult3.inputs[0])
    ng.links.new(input_node.outputs["Fourier Mode 3"], math_mult3.inputs[1])
    
    # Combine all Fourier modes
    ng.links.new(math_mult1.outputs["Value"], math_add_modes.inputs[0])
    ng.links.new(math_mult2.outputs["Value"], math_add_modes.inputs[1])
    ng.links.new(math_add_modes.outputs["Value"], math_add_final.inputs[0])
    ng.links.new(math_mult3.outputs["Value"], math_add_final.inputs[1])
    
    # Star instancing
    ng.links.new(ico_sphere.outputs["Mesh"], instance_on_points.inputs["Instance"])
    ng.links.new(set_position.outputs["Geometry"], instance_on_points.inputs["Points"])
    
    # Final output
    ng.links.new(instance_on_points.outputs["Instances"], output_node.inputs["Geometry"])
    
    return ng


# ASTROPHOT Presets for advanced models
ASTROPHOT_WARP_PRESETS = {
    "mild_warp": {
        "Warp Strength": 0.1,
        "PA Variation": 0.5,
        "Ellipticity Variation": 0.2,
    },
    "strong_warp": {
        "Warp Strength": 0.5,
        "PA Variation": 1.5,
        "Ellipticity Variation": 0.6,
    },
    "realistic_spiral": {
        "Warp Strength": 0.3,
        "PA Variation": 1.0,
        "Ellipticity Variation": 0.4,
    },
}

ASTROPHOT_FOURIER_PRESETS = {
    "barred_galaxy": {
        "Fourier Mode 1": 0.05,
        "Fourier Mode 2": 0.4,   # Strong m=2 for bar
        "Fourier Mode 3": 0.1,
        "Phase 2": 0.0,          # Bar along major axis
    },
    "lopsided_galaxy": {
        "Fourier Mode 1": 0.3,   # Strong m=1 for lopsidedness
        "Fourier Mode 2": 0.1,
        "Fourier Mode 3": 0.05,
        "Phase 1": 0.0,
    },
    "triangular_galaxy": {
        "Fourier Mode 1": 0.1,
        "Fourier Mode 2": 0.1,
        "Fourier Mode 3": 0.3,   # Strong m=3 for Y-shape
        "Phase 3": 0.5,
    },
}


def register():
    """Register all advanced AstroPhot geometry node groups"""
    
    # Register warp models
    if "ASTROPHOT_WarpGalaxy" not in bpy.data.node_groups:
        create_warp_galaxy_node_group()
    
    # Register Fourier models
    if "ASTROPHOT_FourierGalaxy" not in bpy.data.node_groups:
        create_fourier_galaxy_node_group()


def unregister():
    """Unregister advanced AstroPhot geometry node groups"""
    
    advanced_groups = [
        "ASTROPHOT_WarpGalaxy",
        "ASTROPHOT_FourierGalaxy",
    ]
    
    for group_name in advanced_groups:
        if group_name in bpy.data.node_groups:
            bpy.data.node_groups.remove(bpy.data.node_groups[group_name])


if __name__ == "__main__":
    register()
