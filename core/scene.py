"""
Scene Setup for Astronomical Extension
=====================================

Modern scene setup using Blender 4.4 APIs without external config dependencies.
Following clean refactoring principles from graph knowledge.
"""

import bpy
from mathutils import Vector

# ==========================================
# Default Configuration (no external dependencies)
# ==========================================

DEFAULT_CONFIG = {
    "background_color": (0.0, 0.0, 0.0, 1.0),  # Pure black space
    "world_strength": 0.0,  # No ambient light
    "units_system": "METRIC",
    "units_scale": 1.0,
}

PRESET_CONFIGS = {
    "scientific": {
        "background_color": (0.0, 0.0, 0.0, 1.0),
        "world_strength": 0.0,
        "color_temperature": 6500,  # Daylight white balance
    },
    "cinematic": {
        "background_color": (0.02, 0.02, 0.05, 1.0),  # Deep blue-black
        "world_strength": 0.1,
        "color_temperature": 5000,  # Slightly warm
    },
    "publication": {
        "background_color": (1.0, 1.0, 1.0, 1.0),  # White background
        "world_strength": 0.5,
        "color_temperature": 6500,  # Neutral white
    },
}

# ==========================================
# Scene Setup Functions
# ==========================================

def setup_scene(scene=None, preset="scientific", config=None):
    """
    Prepare the Blender scene with modern 4.4 APIs
    
    Args:
        scene: bpy.types.Scene (default: current scene)
        preset: Preset name ('scientific', 'cinematic', 'publication')
        config: Dict to override preset values
    
    Returns:
        dict: Applied configuration
    """
    if scene is None:
        scene = bpy.context.scene
    
    # Get configuration
    final_config = DEFAULT_CONFIG.copy()
    if preset in PRESET_CONFIGS:
        final_config.update(PRESET_CONFIGS[preset])
    if config:
        final_config.update(config)
    
    # Setup units
    scene.unit_settings.system = final_config.get("units_system", "METRIC")
    scene.unit_settings.scale_length = final_config.get("units_scale", 1.0)
    
    # Setup world
    _setup_world(scene, final_config)
    
    # Clean scene (optional - only remove mesh objects)
    if final_config.get("clean_scene", True):
        _clean_scene_objects(scene)
    
    # Setup timeline for animations
    scene.frame_start = 1
    scene.frame_end = 250
    scene.frame_current = 1
    
    print(f"Scene setup completed with preset: {preset}")
    return final_config

def _setup_world(scene, config):
    """Setup world background and environment"""
    # Create or get world
    if scene.world is None:
        scene.world = bpy.data.worlds.new("AstronomicalWorld")
    
    world = scene.world
    world.use_nodes = True
    
    # Get nodes
    nodes = world.node_tree.nodes
    links = world.node_tree.links
    
    # Clear existing nodes
    nodes.clear()
    
    # Create background shader (modern Blender 4.4 approach)
    background_node = nodes.new(type='ShaderNodeBackground')
    background_node.location = (0, 0)
    
    # Set background color
    bg_color = config.get("background_color", (0.0, 0.0, 0.0, 1.0))
    background_node.inputs['Color'].default_value = bg_color
    background_node.inputs['Strength'].default_value = config.get("world_strength", 0.0)
    
    # Create output node
    output_node = nodes.new(type='ShaderNodeOutputWorld')
    output_node.location = (200, 0)
    
    # Link nodes
    links.new(background_node.outputs['Background'], output_node.inputs['Surface'])
    
    # Optional: Add starfield texture for cinematic preset
    if config.get("add_starfield", False):
        _add_starfield_texture(world, nodes, links, background_node)

def _add_starfield_texture(world, nodes, links, background_node):
    """Add procedural starfield texture"""
    # Create noise texture for stars
    noise_node = nodes.new(type='ShaderNodeTexNoise')
    noise_node.location = (-400, 0)
    noise_node.inputs['Scale'].default_value = 1000.0
    noise_node.inputs['Detail'].default_value = 0.0
    noise_node.inputs['Roughness'].default_value = 0.0
    
    # Create ColorRamp for star threshold
    colorramp_node = nodes.new(type='ShaderNodeValToRGB')
    colorramp_node.location = (-200, 0)
    
    # Configure color ramp for sparse stars
    colorramp = colorramp_node.color_ramp
    colorramp.elements[0].position = 0.99  # Very high threshold
    colorramp.elements[0].color = (0, 0, 0, 1)  # Black
    colorramp.elements[1].position = 1.0
    colorramp.elements[1].color = (1, 1, 1, 1)  # White stars
    
    # Mix with background
    mix_node = nodes.new(type='ShaderNodeMix')
    mix_node.location = (-100, 0)
    mix_node.data_type = 'RGBA'
    mix_node.blend_type = 'ADD'
    
    # Link starfield
    links.new(noise_node.outputs['Fac'], colorramp_node.inputs['Fac'])
    links.new(colorramp_node.outputs['Color'], mix_node.inputs['Color2'])
    links.new(mix_node.outputs['Result'], background_node.inputs['Color'])

def _clean_scene_objects(scene, keep_types=None):
    """
    Clean scene objects, keeping only specified types
    
    Args:
        scene: Target scene
        keep_types: Set of object types to keep (default: cameras and lights)
    """
    if keep_types is None:
        keep_types = {"CAMERA", "LIGHT"}
    
    # Get objects to remove
    objects_to_remove = [
        obj for obj in scene.objects 
        if obj.type not in keep_types
    ]
    
    # Remove objects
    for obj in objects_to_remove:
        bpy.data.objects.remove(obj, do_unlink=True)
    
    print(f"Cleaned {len(objects_to_remove)} objects from scene")

def setup_astronomical_view_settings(scene=None):
    """Setup view settings optimized for astronomical rendering"""
    if scene is None:
        scene = bpy.context.scene
    
    # Setup color management for astronomical data
    scene.view_settings.view_transform = 'Standard'
    scene.view_settings.look = 'None'
    scene.view_settings.exposure = 0.0
    scene.view_settings.gamma = 1.0
    
    # Setup sequencer color space if available
    if hasattr(scene.sequencer_colorspace_settings, 'name'):
        scene.sequencer_colorspace_settings.name = 'sRGB'

def reset_scene_to_defaults():
    """Reset scene to factory defaults for astronomical work"""
    # Use factory settings (modern Blender 4.4 approach)
    bpy.ops.wm.read_factory_settings(use_empty=True)
    
    # Apply astronomical setup
    setup_scene(preset="scientific")
    
    print("Scene reset to astronomical defaults")

# ==========================================
# Utility Functions
# ==========================================

def get_scene_info(scene=None):
    """Get information about current scene setup"""
    if scene is None:
        scene = bpy.context.scene
    
    return {
        "name": scene.name,
        "frame_range": (scene.frame_start, scene.frame_end),
        "current_frame": scene.frame_current,
        "render_engine": scene.render.engine,
        "units": scene.unit_settings.system,
        "objects_count": len(scene.objects),
        "world": scene.world.name if scene.world else None,
    }

def validate_scene_for_astronomy():
    """Validate that scene is properly configured for astronomical rendering"""
    scene = bpy.context.scene
    issues = []
    
    # Check render engine
    if scene.render.engine not in ['CYCLES', 'BLENDER_EEVEE_NEXT']:
        issues.append(f"Render engine '{scene.render.engine}' not optimal for astronomy")
    
    # Check world setup
    if not scene.world:
        issues.append("No world material assigned")
    elif not scene.world.use_nodes:
        issues.append("World material not using nodes")
    
    # Check units
    if scene.unit_settings.system != 'METRIC':
        issues.append("Units not set to METRIC")
    
    return {
        "valid": len(issues) == 0,
        "issues": issues,
        "scene_info": get_scene_info(scene)
    }
