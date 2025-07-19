"""
Rendering Configuration for Astronomical Extension
=================================================

Modern rendering setup using Blender 4.4 EEVEE Next and Cycles.
Following clean refactoring principles - no external config dependencies.
"""

import bpy
import os
from pathlib import Path
from typing import Optional, Tuple, Dict, Any

# ==========================================
# Default Rendering Configuration
# ==========================================

DEFAULT_RENDER_CONFIG = {
    "engine": "BLENDER_EEVEE_NEXT",  # Modern Blender 4.4 default
    "resolution": (1920, 1080),
    "samples": 64,
    "file_format": "PNG",
    "output_path": "astronomical_render.png",
    "color_management": {
        "display_device": "sRGB",
        "view_transform": "Standard",  # Modern default
        "look": "None",
    }
}

PRESET_CONFIGS = {
    "preview": {
        "engine": "BLENDER_EEVEE_NEXT",
        "resolution": (640, 480),
        "samples": 16,
        "file_format": "JPEG",
    },
    "final": {
        "engine": "CYCLES", 
        "resolution": (3840, 2160),  # 4K
        "samples": 256,
        "file_format": "OPEN_EXR",
    },
    "animation": {
        "engine": "BLENDER_EEVEE_NEXT",
        "resolution": (1920, 1080),
        "samples": 32,
        "file_format": "FFMPEG",
    }
}

# ==========================================
# Rendering Setup Functions
# ==========================================

def setup_rendering(scene=None, preset="default", config=None):
    """
    Configure rendering settings for astronomical visualization
    
    Args:
        scene: Target scene (default: current scene)
        preset: Preset name ('preview', 'final', 'animation', 'default')
        config: Dict to override preset values
    
    Returns:
        dict: Applied configuration
    """
    if scene is None:
        scene = bpy.context.scene
    
    # Get configuration
    final_config = DEFAULT_RENDER_CONFIG.copy()
    if preset in PRESET_CONFIGS:
        final_config.update(PRESET_CONFIGS[preset])
    if config:
        final_config.update(config)
    
    # Apply engine settings
    _setup_render_engine(scene, final_config)
    
    # Apply resolution and format
    _setup_output_settings(scene, final_config)
    
    # Apply color management
    _setup_color_management(scene, final_config)
    
    # Apply engine-specific settings
    _setup_engine_specific_settings(scene, final_config)
    
    print(f"Rendering configured with engine: {scene.render.engine}")
    return final_config

def _setup_render_engine(scene, config):
    """Setup render engine and basic settings"""
    engine = config.get("engine", "BLENDER_EEVEE_NEXT")
    
    # Validate engine for Blender 4.4
    if engine not in ['CYCLES', 'BLENDER_EEVEE_NEXT', 'WORKBENCH']:
        print(f"Warning: Unknown engine '{engine}', using EEVEE Next")
        engine = 'BLENDER_EEVEE_NEXT'
    
    scene.render.engine = engine
    
    # Basic render settings
    scene.render.resolution_percentage = 100
    scene.render.use_persistent_data = True
    
    # Frame settings
    if not hasattr(scene, '_astro_frame_configured'):
        scene.frame_start = 1
        scene.frame_end = 1
        scene.frame_current = 1
        scene._astro_frame_configured = True

def _setup_output_settings(scene, config):
    """Setup output resolution and file format"""
    # Resolution
    resolution = config.get("resolution", (1920, 1080))
    scene.render.resolution_x = resolution[0]
    scene.render.resolution_y = resolution[1]
    
    # File format
    file_format = config.get("file_format", "PNG")
    scene.render.image_settings.file_format = file_format
    
    # Format-specific settings
    if file_format == "PNG":
        scene.render.image_settings.color_mode = 'RGBA'
        scene.render.image_settings.color_depth = '16'
        scene.render.image_settings.compression = 90
    elif file_format == "OPEN_EXR":
        scene.render.image_settings.color_mode = 'RGBA'
        scene.render.image_settings.color_depth = '32'
        scene.render.image_settings.exr_codec = 'ZIP'
    elif file_format == "JPEG":
        scene.render.image_settings.color_mode = 'RGB'
        scene.render.image_settings.quality = 90
    
    # Output path
    output_path = config.get("output_path", "astronomical_render.png")
    scene.render.filepath = _ensure_output_directory(output_path)

def _setup_color_management(scene, config):
    """Setup color management for astronomical rendering"""
    color_config = config.get("color_management", {})
    
    # Display settings
    scene.display_settings.display_device = color_config.get("display_device", "sRGB")
    
    # View settings
    scene.view_settings.view_transform = color_config.get("view_transform", "Standard")
    scene.view_settings.look = color_config.get("look", "None")
    scene.view_settings.exposure = color_config.get("exposure", 0.0)
    scene.view_settings.gamma = color_config.get("gamma", 1.0)

def _setup_engine_specific_settings(scene, config):
    """Setup engine-specific rendering settings"""
    samples = config.get("samples", 64)
    
    if scene.render.engine == "CYCLES":
        _setup_cycles_settings(scene, config, samples)
    elif scene.render.engine == "BLENDER_EEVEE_NEXT":
        _setup_eevee_next_settings(scene, config, samples)

def _setup_cycles_settings(scene, config, samples):
    """Setup Cycles render engine settings"""
    cycles = scene.cycles
    
    # Sampling
    cycles.samples = samples
    cycles.preview_samples = max(16, samples // 4)
    
    # Performance
    cycles.use_persistent_data = True
    
    # Feature set
    cycles.feature_set = 'SUPPORTED'
    
    # Device (prefer GPU if available)
    preferences = bpy.context.preferences
    cycles_prefs = preferences.addons['cycles'].preferences
    if cycles_prefs.compute_device_type != 'NONE':
        cycles.device = 'GPU'
    else:
        cycles.device = 'CPU'
    
    # Light paths (optimized for astronomy)
    cycles.max_bounces = 4
    cycles.diffuse_bounces = 2
    cycles.glossy_bounces = 2
    cycles.transmission_bounces = 2
    cycles.volume_bounces = 0  # Usually not needed for astronomy
    
    # Denoising
    cycles.use_denoising = True
    cycles.denoiser = 'OPENIMAGEDENOISE'

def _setup_eevee_next_settings(scene, config, samples):
    """Setup EEVEE Next render engine settings (Blender 4.4)"""
    eevee = scene.eevee
    
    # Sampling (modern EEVEE Next API)
    if hasattr(eevee, 'taa_render_samples'):
        eevee.taa_render_samples = samples
    
    # Quality settings
    if hasattr(eevee, 'use_ssr'):
        eevee.use_ssr = True  # Screen space reflections
        eevee.ssr_quality = 0.25  # Balance quality/performance
    
    if hasattr(eevee, 'use_volumetric_lights'):
        eevee.use_volumetric_lights = True
        eevee.volumetric_tile_size = '2'  # Good for astronomical effects
    
    # Motion blur (for animations)
    if hasattr(eevee, 'use_motion_blur'):
        eevee.use_motion_blur = False  # Usually off for astronomy
    
    # Ambient occlusion
    if hasattr(eevee, 'use_gtao'):
        eevee.use_gtao = True
        eevee.gtao_distance = 1.0

def _ensure_output_directory(output_path):
    """Ensure output directory exists and return absolute path"""
    path = Path(output_path)
    
    # Make absolute if relative
    if not path.is_absolute():
        # Use project directory or current working directory
        base_dir = Path.cwd()
        if hasattr(bpy.context, 'blend_data') and bpy.data.filepath:
            base_dir = Path(bpy.data.filepath).parent
        path = base_dir / path
    
    # Create directory
    path.parent.mkdir(parents=True, exist_ok=True)
    
    return str(path)

# ==========================================
# Main Rendering Function
# ==========================================

def render_astronomical_scene(
    output_path: Optional[str] = None,
    render_engine: Optional[str] = None, 
    resolution: Optional[Tuple[int, int]] = None,
    samples: Optional[int] = None,
    preset: str = "default",
    **kwargs
) -> Dict[str, Any]:
    """
    Render astronomical scene with professional settings
    
    Args:
        output_path: Output file path
        render_engine: Render engine ('CYCLES', 'BLENDER_EEVEE_NEXT')
        resolution: Render resolution (width, height)
        samples: Number of samples
        preset: Configuration preset
        **kwargs: Additional render parameters
    
    Returns:
        dict: Render information and results
    """
    scene = bpy.context.scene
    
    # Build configuration
    config = {}
    if render_engine:
        config["engine"] = render_engine
    if resolution:
        config["resolution"] = resolution
    if samples:
        config["samples"] = samples
    if output_path:
        config["output_path"] = output_path
    
    # Apply additional parameters
    config.update(kwargs)
    
    # Setup rendering
    final_config = setup_rendering(scene, preset, config)
    
    # Pre-render checks
    render_info = _validate_scene_for_render(scene)
    if not render_info["valid"]:
        print("Warning: Scene validation issues found:")
        for issue in render_info["issues"]:
            print(f"  - {issue}")
    
    # Perform render
    start_time = time.time()
    
    try:
        print(f"Starting render with {scene.render.engine}...")
        print(f"Resolution: {scene.render.resolution_x}x{scene.render.resolution_y}")
        print(f"Output: {scene.render.filepath}")
        
        # Render single frame
        bpy.ops.render.render(write_still=True)
        
        render_time = time.time() - start_time
        
        result = {
            "success": True,
            "output_path": scene.render.filepath,
            "render_time": render_time,
            "engine": scene.render.engine,
            "resolution": (scene.render.resolution_x, scene.render.resolution_y),
            "config": final_config,
        }
        
        print(f"✅ Render completed in {render_time:.2f} seconds")
        print(f"Output saved to: {result['output_path']}")
        
        return result
        
    except Exception as e:
        print(f"❌ Render failed: {e}")
        return {
            "success": False,
            "error": str(e),
            "config": final_config,
        }

def _validate_scene_for_render(scene):
    """Validate scene is ready for astronomical rendering"""
    issues = []
    
    # Check for objects
    if len(scene.objects) == 0:
        issues.append("Scene contains no objects")
    
    # Check for camera
    if not scene.camera:
        issues.append("No active camera in scene")
    
    # Check for lights (astronomy may intentionally have no lights)
    lights = [obj for obj in scene.objects if obj.type == 'LIGHT']
    if len(lights) == 0:
        issues.append("No lights in scene (may be intentional for astronomy)")
    
    # Check render engine availability
    if scene.render.engine not in ['CYCLES', 'BLENDER_EEVEE_NEXT', 'WORKBENCH']:
        issues.append(f"Unknown render engine: {scene.render.engine}")
    
    return {
        "valid": len(issues) == 0,
        "issues": issues,
        "object_count": len(scene.objects),
        "light_count": len(lights),
        "has_camera": scene.camera is not None,
    }

# ==========================================
# Utility Functions
# ==========================================

import time

def get_render_info(scene=None):
    """Get current render configuration info"""
    if scene is None:
        scene = bpy.context.scene
    
    return {
        "engine": scene.render.engine,
        "resolution": (scene.render.resolution_x, scene.render.resolution_y),
        "file_format": scene.render.image_settings.file_format,
        "output_path": scene.render.filepath,
        "samples": _get_current_samples(scene),
        "color_space": scene.view_settings.view_transform,
    }

def _get_current_samples(scene):
    """Get current sample count for active render engine"""
    if scene.render.engine == "CYCLES":
        return scene.cycles.samples
    elif scene.render.engine == "BLENDER_EEVEE_NEXT":
        if hasattr(scene.eevee, 'taa_render_samples'):
            return scene.eevee.taa_render_samples
    return "Unknown"

def optimize_for_headless():
    """Optimize settings for headless rendering"""
    # Disable unnecessary features for headless operation
    preferences = bpy.context.preferences
    
    # Disable add-ons that might cause issues
    preferences.view.show_splash = False
    
    # Use CPU rendering for better headless compatibility
    scene = bpy.context.scene
    if scene.render.engine == "CYCLES":
        scene.cycles.device = 'CPU'
    
    print("Render settings optimized for headless operation")

def setup_astronomical_post_processing(scene=None):
    """Setup post-processing nodes for astronomical rendering"""
    if scene is None:
        scene = bpy.context.scene
    
    # Enable compositor
    scene.use_nodes = True
    
    # Clear existing nodes
    scene.node_tree.nodes.clear()
    
    # Create render layer input
    render_layers = scene.node_tree.nodes.new('CompositorNodeRLayers')
    render_layers.location = (0, 0)
    
    # Create output
    composite = scene.node_tree.nodes.new('CompositorNodeComposite')
    composite.location = (300, 0)
    
    # Link basic setup
    scene.node_tree.links.new(
        render_layers.outputs['Image'],
        composite.inputs['Image']
    )
    
    # Optional: Add astronomical-specific post-processing
    # This could include star glow, color grading, etc.
    
    print("Astronomical post-processing setup completed")
