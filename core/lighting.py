"""
Lighting Setup for Astronomical Extension
========================================

Modern lighting configuration for astronomical visualization.
No external config dependencies - following clean refactoring principles.
"""

import bpy
from mathutils import Vector
from typing import Optional, List, Dict, Any

# ==========================================
# Default Lighting Configuration
# ==========================================

DEFAULT_LIGHTING_CONFIG = {
    "preset": "space",
    "remove_existing": True,
    "ambient_strength": 0.01,
}

LIGHTING_PRESETS = {
    "space": {
        "description": "Minimal lighting for space scenes",
        "lights": [
            {
                "name": "Starlight",
                "type": "SUN",
                "energy": 0.05,
                "color": (0.9, 0.95, 1.0),  # Cool blue-white
                "location": (10, -10, 10),
                "rotation": (0.785, 0, 0.785),  # 45° angles
            }
        ]
    },
    "galaxy": {
        "description": "Galactic center glow with ambient starlight",
        "lights": [
            {
                "name": "Galactic_Glow", 
                "type": "POINT",
                "energy": 100,
                "color": (1.0, 0.8, 0.6),  # Warm galactic center
                "location": (0, 0, 0),
                "falloff_type": "INVERSE_SQUARE",
            },
            {
                "name": "Ambient_Starlight",
                "type": "SUN", 
                "energy": 0.01,
                "color": (0.9, 0.95, 1.0),
                "location": (50, -50, 50),
            }
        ]
    },
    "planetary": {
        "description": "Solar system lighting with main star",
        "lights": [
            {
                "name": "Central_Star",
                "type": "POINT",
                "energy": 1000,
                "color": (1.0, 1.0, 0.9),  # Solar white
                "location": (0, 0, 0),
                "shadow_soft_size": 0.5,
            },
            {
                "name": "Ambient_Space",
                "type": "SUN",
                "energy": 0.005,
                "color": (0.8, 0.9, 1.0),
                "location": (100, -100, 100),
            }
        ]
    },
    "nebula": {
        "description": "Colorful nebula lighting with emission",
        "lights": [
            {
                "name": "Nebula_Core",
                "type": "AREA",
                "energy": 50,
                "color": (1.0, 0.3, 0.8),  # Magenta
                "location": (0, 0, 5),
                "size": 10,
                "shape": "DISK",
            },
            {
                "name": "Blue_Region",
                "type": "AREA", 
                "energy": 30,
                "color": (0.2, 0.5, 1.0),  # Blue
                "location": (-10, 10, 0),
                "size": 15,
                "shape": "DISK",
            }
        ]
    },
    "scientific": {
        "description": "Neutral lighting for scientific visualization",
        "lights": [
            {
                "name": "Key_Light",
                "type": "SUN",
                "energy": 2.0,
                "color": (1.0, 1.0, 1.0),  # Pure white
                "location": (10, -10, 10),
            },
            {
                "name": "Fill_Light",
                "type": "SUN",
                "energy": 0.5,
                "color": (0.9, 0.95, 1.0),  # Slightly cool
                "location": (-5, 5, 8),
            }
        ]
    }
}

# ==========================================
# Lighting Setup Functions
# ==========================================

def setup_lighting(scene=None, preset="space", config=None):
    """
    Setup lighting for astronomical visualization
    
    Args:
        scene: Target scene (default: current scene)
        preset: Lighting preset name
        config: Dict to override preset values
    
    Returns:
        list: Created light objects
    """
    if scene is None:
        scene = bpy.context.scene
    
    # Get configuration
    final_config = DEFAULT_LIGHTING_CONFIG.copy()
    if config:
        final_config.update(config)
    
    # Remove existing lights if requested
    if final_config.get("remove_existing", True):
        _remove_existing_lights(scene)
    
    # Create lights from preset
    if preset not in LIGHTING_PRESETS:
        print(f"Warning: Unknown lighting preset '{preset}', using 'space'")
        preset = "space"
    
    lights = _create_lights_from_preset(scene, preset)
    
    # Setup world ambient if specified
    _setup_world_ambient(scene, final_config)
    
    print(f"Lighting setup completed: {preset} ({len(lights)} lights)")
    return lights

def _remove_existing_lights(scene):
    """Remove all existing lights from scene"""
    lights_to_remove = [obj for obj in scene.objects if obj.type == "LIGHT"]
    
    for light_obj in lights_to_remove:
        bpy.data.objects.remove(light_obj, do_unlink=True)
    
    if lights_to_remove:
        print(f"Removed {len(lights_to_remove)} existing lights")

def _create_lights_from_preset(scene, preset_name):
    """Create lights from preset configuration"""
    preset = LIGHTING_PRESETS[preset_name]
    created_lights = []
    
    for light_config in preset["lights"]:
        light_obj = _create_single_light(scene, light_config)
        created_lights.append(light_obj)
    
    return created_lights

def _create_single_light(scene, light_config):
    """Create a single light from configuration"""
    # Create light data
    light_data = bpy.data.lights.new(
        name=light_config["name"],
        type=light_config["type"]
    )
    
    # Create light object
    light_obj = bpy.data.objects.new(light_config["name"], light_data)
    scene.collection.objects.link(light_obj)
    
    # Configure light properties
    _configure_light_properties(light_data, light_config)
    
    # Position light
    light_obj.location = Vector(light_config.get("location", (0, 0, 10)))
    
    # Set rotation if specified
    if "rotation" in light_config:
        light_obj.rotation_euler = light_config["rotation"]
    
    return light_obj

def _configure_light_properties(light_data, config):
    """Configure light data properties"""
    # Basic properties
    light_data.energy = config.get("energy", 1.0)
    light_data.color = config.get("color", (1.0, 1.0, 1.0))
    
    # Type-specific properties
    if light_data.type == "SUN":
        light_data.angle = config.get("angle", 0.00918)  # Sun's angular size
        
    elif light_data.type == "POINT":
        light_data.shadow_soft_size = config.get("shadow_soft_size", 0.1)
        
    elif light_data.type == "SPOT":
        light_data.spot_size = config.get("spot_size", 0.785)  # 45°
        light_data.spot_blend = config.get("spot_blend", 0.15)
        
    elif light_data.type == "AREA":
        light_data.size = config.get("size", 1.0)
        if "shape" in config:
            light_data.shape = config["shape"]
        if light_data.shape == "RECTANGLE":
            light_data.size_y = config.get("size_y", light_data.size)

def _setup_world_ambient(scene, config):
    """Setup world ambient lighting"""
    ambient_strength = config.get("ambient_strength", 0.01)
    
    if ambient_strength > 0 and scene.world:
        # Ensure world uses nodes
        scene.world.use_nodes = True
        
        # Find background node
        nodes = scene.world.node_tree.nodes
        bg_node = None
        for node in nodes:
            if node.type == 'BACKGROUND':
                bg_node = node
                break
        
        if bg_node:
            # Add minimal ambient to background
            current_strength = bg_node.inputs['Strength'].default_value
            bg_node.inputs['Strength'].default_value = max(current_strength, ambient_strength)

# ==========================================
# Specialized Lighting Functions  
# ==========================================

def create_star_lighting(star_temperature=5778, luminosity=1.0, location=(0, 0, 0)):
    """
    Create realistic star lighting based on temperature
    
    Args:
        star_temperature: Star temperature in Kelvin (default: Sun = 5778K)
        luminosity: Relative luminosity (1.0 = Sun)
        location: Star position
    
    Returns:
        bpy.types.Object: Star light object
    """
    scene = bpy.context.scene
    
    # Calculate color from temperature (simplified blackbody)
    color = _temperature_to_rgb(star_temperature)
    
    # Create point light
    light_data = bpy.data.lights.new("Star_Light", "POINT")
    light_data.energy = 1000 * luminosity
    light_data.color = color
    light_data.shadow_soft_size = 0.01  # Sharp shadows for distant star
    
    # Create object
    light_obj = bpy.data.objects.new("Star_Light", light_data)
    light_obj.location = Vector(location)
    scene.collection.objects.link(light_obj)
    
    return light_obj

def _temperature_to_rgb(temperature):
    """Convert color temperature to RGB (simplified)"""
    # Simplified blackbody color conversion
    if temperature < 3500:
        return (1.0, 0.4, 0.2)  # Red dwarf
    elif temperature < 5000:
        return (1.0, 0.8, 0.6)  # Orange
    elif temperature < 6000:
        return (1.0, 1.0, 0.9)  # Yellow-white (Sun)
    elif temperature < 8000:
        return (0.9, 0.95, 1.0)  # White
    else:
        return (0.7, 0.8, 1.0)  # Blue giant

def create_nebula_lighting(center=(0, 0, 0), size=20, color_scheme="emission"):
    """
    Create multi-colored nebula lighting
    
    Args:
        center: Nebula center position
        size: Nebula size
        color_scheme: Color scheme ('emission', 'reflection', 'planetary')
    
    Returns:
        list: Created light objects
    """
    scene = bpy.context.scene
    lights = []
    
    color_schemes = {
        "emission": [(1.0, 0.2, 0.3), (0.2, 0.8, 0.3), (0.3, 0.2, 1.0)],  # H-alpha, OIII, other
        "reflection": [(0.7, 0.8, 1.0), (0.9, 0.9, 1.0)],  # Blue reflection
        "planetary": [(0.2, 1.0, 0.8), (1.0, 0.8, 0.2)],  # Teal and orange
    }
    
    colors = color_schemes.get(color_scheme, color_schemes["emission"])
    
    # Create multiple area lights with different colors
    for i, color in enumerate(colors):
        light_data = bpy.data.lights.new(f"Nebula_Region_{i}", "AREA")
        light_data.energy = 20
        light_data.color = color
        light_data.size = size
        light_data.shape = "DISK"
        
        # Position lights around center
        import math
        angle = (i / len(colors)) * 2 * math.pi
        offset_x = math.cos(angle) * size * 0.3
        offset_y = math.sin(angle) * size * 0.3
        
        light_obj = bpy.data.objects.new(f"Nebula_Region_{i}", light_data)
        light_obj.location = Vector((center[0] + offset_x, center[1] + offset_y, center[2]))
        scene.collection.objects.link(light_obj)
        
        lights.append(light_obj)
    
    return lights

# ==========================================
# Utility Functions
# ==========================================

def get_lighting_info(scene=None):
    """Get information about current lighting setup"""
    if scene is None:
        scene = bpy.context.scene
    
    lights = [obj for obj in scene.objects if obj.type == "LIGHT"]
    
    light_info = []
    for light_obj in lights:
        light_data = light_obj.data
        light_info.append({
            "name": light_obj.name,
            "type": light_data.type,
            "energy": light_data.energy,
            "color": tuple(light_data.color),
            "location": tuple(light_obj.location),
        })
    
    return {
        "light_count": len(lights),
        "lights": light_info,
        "world_background": scene.world.node_tree.nodes["Background"].inputs["Strength"].default_value if scene.world and scene.world.use_nodes else 0,
    }

def optimize_lighting_for_engine(scene=None, engine="BLENDER_EEVEE_NEXT"):
    """Optimize lighting settings for specific render engine"""
    if scene is None:
        scene = bpy.context.scene
    
    lights = [obj for obj in scene.objects if obj.type == "LIGHT"]
    
    if engine == "BLENDER_EEVEE_NEXT":
        # EEVEE Next optimizations
        for light_obj in lights:
            light_data = light_obj.data
            # Enable contact shadows for realism
            if hasattr(light_data, 'use_contact_shadow'):
                light_data.use_contact_shadow = True
            
    elif engine == "CYCLES":
        # Cycles optimizations
        for light_obj in lights:
            light_data = light_obj.data
            # Increase samples for area lights
            if light_data.type == "AREA":
                light_data.cycles.samples = 8
    
    print(f"Lighting optimized for {engine}")

def create_lighting_preset_from_scene(preset_name):
    """Create new lighting preset from current scene"""
    scene = bpy.context.scene
    lights = [obj for obj in scene.objects if obj.type == "LIGHT"]
    
    preset_config = {
        "description": f"Custom preset: {preset_name}",
        "lights": []
    }
    
    for light_obj in lights:
        light_data = light_obj.data
        light_config = {
            "name": light_obj.name,
            "type": light_data.type,
            "energy": light_data.energy,
            "color": tuple(light_data.color),
            "location": tuple(light_obj.location),
            "rotation": tuple(light_obj.rotation_euler),
        }
        
        # Add type-specific properties
        if light_data.type == "AREA":
            light_config["size"] = light_data.size
            light_config["shape"] = light_data.shape
        elif light_data.type == "SPOT":
            light_config["spot_size"] = light_data.spot_size
            light_config["spot_blend"] = light_data.spot_blend
        
        preset_config["lights"].append(light_config)
    
    print(f"Created preset '{preset_name}' with {len(lights)} lights")
    return preset_config
