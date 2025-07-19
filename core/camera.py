"""
Camera Setup for Astronomical Extension
======================================

Modern camera configuration for astronomical visualization.
No external config dependencies - following clean refactoring principles.
"""

import bpy
from mathutils import Vector, Matrix
from typing import Optional, Tuple, Dict, Any

# ==========================================
# Default Camera Configuration
# ==========================================

DEFAULT_CAMERA_CONFIG = {
    "location": (0, -15, 5),
    "look_at": (0, 0, 0), 
    "focal_length": 50,
    "sensor_width": 36,  # Full frame sensor
    "clip_start": 0.1,
    "clip_end": 10000,  # Large for cosmic scales
}

PRESET_CONFIGS = {
    "stellar_field": {
        "location": (0, -50, 20),
        "look_at": (0, 0, 0),
        "focal_length": 85,  # Telephoto for stars
        "clip_end": 10000,
    },
    "galaxy": {
        "location": (50, -50, 30),
        "look_at": (0, 0, 0),
        "focal_length": 35,  # Wide angle
        "clip_end": 1000,
    },
    "planetary": {
        "location": (0, -30, 10),
        "look_at": (0, 0, 0),
        "focal_length": 135,  # Telephoto for planets
        "clip_end": 500,
    },
    "cosmic_web": {
        "location": (200, -200, 100),
        "look_at": (0, 0, 0),
        "focal_length": 24,  # Ultra wide
        "clip_end": 50000,
    }
}

# ==========================================
# Camera Setup Functions
# ==========================================

def setup_camera(scene=None, preset="default", config=None):
    """
    Setup camera for astronomical visualization
    
    Args:
        scene: Target scene (default: current scene)
        preset: Camera preset ('stellar_field', 'galaxy', 'planetary', 'cosmic_web')
        config: Dict to override preset values
    
    Returns:
        bpy.types.Object: Camera object
    """
    if scene is None:
        scene = bpy.context.scene
    
    # Get configuration
    final_config = DEFAULT_CAMERA_CONFIG.copy()
    if preset in PRESET_CONFIGS:
        final_config.update(PRESET_CONFIGS[preset])
    if config:
        final_config.update(config)
    
    # Find or create camera
    camera_obj = _get_or_create_camera(scene)
    
    # Configure camera
    _configure_camera_settings(camera_obj, final_config)
    _position_camera(camera_obj, final_config)
    
    # Set as active camera
    scene.camera = camera_obj
    
    print(f"Camera configured for {preset} view")
    return camera_obj

def _get_or_create_camera(scene):
    """Get existing camera or create new one"""
    # Look for existing camera
    for obj in scene.objects:
        if obj.type == "CAMERA":
            return obj
    
    # Create new camera
    camera_data = bpy.data.cameras.new("AstronomicalCamera")
    camera_obj = bpy.data.objects.new("AstronomicalCamera", camera_data)
    scene.collection.objects.link(camera_obj)
    
    return camera_obj

def _configure_camera_settings(camera_obj, config):
    """Configure camera data settings"""
    camera_data = camera_obj.data
    
    # Lens settings
    camera_data.lens = config.get("focal_length", 50)
    camera_data.sensor_width = config.get("sensor_width", 36)
    
    # Clipping
    camera_data.clip_start = config.get("clip_start", 0.1)
    camera_data.clip_end = config.get("clip_end", 10000)
    
    # Type (perspective for astronomy)
    camera_data.type = 'PERSP'
    
    # Depth of field (usually off for astronomy)
    camera_data.dof.use_dof = config.get("use_dof", False)
    if camera_data.dof.use_dof:
        camera_data.dof.aperture_fstop = config.get("f_stop", 2.8)

def _position_camera(camera_obj, config):
    """Position and orient camera"""
    # Set location
    location = config.get("location", (0, -15, 5))
    camera_obj.location = Vector(location)
    
    # Point camera at target
    look_at = Vector(config.get("look_at", (0, 0, 0)))
    _point_camera_at(camera_obj, look_at)

def _point_camera_at(camera_obj, target_location):
    """Point camera at specific location"""
    # Calculate direction vector
    direction = target_location - camera_obj.location
    
    # Create rotation quaternion for -Z axis (camera forward)
    rot_quat = direction.to_track_quat('-Z', 'Y')
    
    # Apply rotation
    camera_obj.rotation_euler = rot_quat.to_euler()

# ==========================================
# Advanced Camera Functions
# ==========================================

def create_orbital_camera(center=(0, 0, 0), radius=20, height=5, orbit_speed=1.0):
    """
    Create camera that orbits around a point (useful for animations)
    
    Args:
        center: Point to orbit around
        radius: Orbital radius
        height: Height above orbital plane
        orbit_speed: Animation speed multiplier
    
    Returns:
        bpy.types.Object: Camera with orbital animation
    """
    scene = bpy.context.scene
    
    # Create camera
    camera_obj = setup_camera(scene, config={
        "location": (radius, 0, height),
        "look_at": center,
        "focal_length": 50
    })
    
    # Add orbital animation
    _add_orbital_animation(camera_obj, center, radius, height, orbit_speed)
    
    return camera_obj

def _add_orbital_animation(camera_obj, center, radius, height, speed):
    """Add orbital animation to camera"""
    import math
    
    # Clear existing animation
    camera_obj.animation_data_clear()
    
    # Create keyframes for orbital motion
    frame_count = 250  # 10 seconds at 25fps
    
    for frame in range(0, frame_count + 1):
        scene = bpy.context.scene
        scene.frame_set(frame)
        
        # Calculate orbital position
        angle = (frame / frame_count) * 2 * math.pi * speed
        
        x = center[0] + radius * math.cos(angle)
        y = center[1] + radius * math.sin(angle)
        z = center[2] + height
        
        # Set position
        camera_obj.location = (x, y, z)
        
        # Point at center
        _point_camera_at(camera_obj, Vector(center))
        
        # Insert keyframes
        camera_obj.keyframe_insert(data_path="location", frame=frame)
        camera_obj.keyframe_insert(data_path="rotation_euler", frame=frame)
    
    # Set interpolation to linear for smooth motion
    if camera_obj.animation_data and camera_obj.animation_data.action:
        for fcurve in camera_obj.animation_data.action.fcurves:
            for keyframe in fcurve.keyframe_points:
                keyframe.interpolation = 'LINEAR'

def get_camera_info(camera_obj=None):
    """Get information about camera configuration"""
    if camera_obj is None:
        scene = bpy.context.scene
        camera_obj = scene.camera
    
    if not camera_obj or camera_obj.type != 'CAMERA':
        return {"error": "No valid camera found"}
    
    camera_data = camera_obj.data
    
    return {
        "name": camera_obj.name,
        "location": tuple(camera_obj.location),
        "rotation": tuple(camera_obj.rotation_euler),
        "focal_length": camera_data.lens,
        "sensor_width": camera_data.sensor_width,
        "clip_start": camera_data.clip_start,
        "clip_end": camera_data.clip_end,
        "type": camera_data.type,
        "dof_enabled": camera_data.dof.use_dof,
    }

def calculate_fov(focal_length, sensor_width=36):
    """
    Calculate field of view from focal length and sensor width
    
    Args:
        focal_length: Lens focal length in mm
        sensor_width: Sensor width in mm (default: 36mm full frame)
    
    Returns:
        float: Field of view in degrees
    """
    import math
    fov_radians = 2 * math.atan(sensor_width / (2 * focal_length))
    return math.degrees(fov_radians)
