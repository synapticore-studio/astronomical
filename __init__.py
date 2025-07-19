"""
Astronomical Blender Extension
=============================

Professional astronomical visualization tools for Blender 4.4+
Clean implementation following strict refactoring principles.

Features:
- Stellar field generation
- Galaxy morphology visualization  
- Cosmic web simulation
- Custom astronomical node groups
- Professional astronomical materials
- Scientific data integration
"""

import bpy
import sys
from pathlib import Path

# Extension Info (for legacy compatibility)
bl_info = {
    "name": "Astronomical Blender Extension",
    "author": "Astronomical Developer",
    "version": (2, 0, 0),
    "blender": (4, 4, 0),
    "location": "View3D > Sidebar > Astronomical",
    "description": "Professional astronomical visualization tools",
    "category": "Add Mesh",
}

# ==========================================
# Vendor Path Setup (following vendorizing approach)
# ==========================================

def setup_vendor_path():
    """Add vendor directory to Python path for dependencies"""
    vendor_path = Path(__file__).parent / "vendor"
    if vendor_path.exists() and str(vendor_path) not in sys.path:
        sys.path.insert(0, str(vendor_path))
        return True
    return False

# Setup vendor imports before any dependency imports (following clean import standards)
setup_vendor_path()

# ==========================================
# Core Imports (all at top level - following graph standards)
# ==========================================

# Core functionality
from .core import (
    setup_camera,
    setup_lighting, 
    setup_rendering,
    setup_scene,
    render_astronomical_scene,
)

# Node system imports (modernized for Blender 4.4)
from .nodes.geometry import register as register_geometry_nodes
from .nodes.geometry import unregister as unregister_geometry_nodes
from .nodes.shader import register as register_shader_nodes  
from .nodes.shader import unregister as unregister_shader_nodes
from .nodes.compositing import register as register_compositing_nodes
from .nodes.compositing import unregister as unregister_compositing_nodes

# Operator imports
from .operators import register as register_operators
from .operators import unregister as unregister_operators

# Utility imports
from .utilities.astronomical_data import (
    STELLAR_CLASSIFICATION,
    GALAXY_TYPES,
    HR_DIAGRAM_PARAMS,
    get_stellar_data,
    get_galaxy_config,
    create_sample_stellar_data,
    validate_hr_diagram_data,
)

# ==========================================
# Extension Registration (modernized for Blender 4.4)
# ==========================================

def register():
    """Register all extension components with modern Blender 4.4 API"""
    try:
        print("Registering Astronomical Blender Extension...")
        
        # Check if running as extension (Blender 4.4 feature)
        if hasattr(bpy.app, 'module') and bpy.app.module:
            print("Running as Blender module - adapting behavior")
        
        # Enable VFX libraries if available (Blender 4.4 feature)
        if hasattr(bpy.utils, 'expose_bundled_modules'):
            try:
                bpy.utils.expose_bundled_modules()
                print("VFX libraries exposed successfully")
            except Exception as e:
                print(f"Could not expose VFX libraries: {e}")
        
        # 1. Register Node Groups (using modern Interface API)
        register_geometry_nodes()
        register_shader_nodes()
        register_compositing_nodes()
        
        # 2. Register Operators
        register_operators()
        
        print("✅ Astronomical Extension registered successfully")
        
    except Exception as e:
        print(f"❌ Extension registration failed: {e}")
        raise

def unregister():
    """Unregister all extension components"""
    try:
        print("Unregistering Astronomical Extension...")
        
        # Unregister in reverse order
        unregister_operators()
        unregister_compositing_nodes()
        unregister_shader_nodes()
        unregister_geometry_nodes()
        
        print("✅ Astronomical Extension unregistered successfully")
        
    except Exception as e:
        print(f"❌ Extension unregistration failed: {e}")

# ==========================================
# Main API Functions
# ==========================================

def create_stellar_field(star_count=1000, distribution="random", scale=1.0, **kwargs):
    """
    Create a stellar field using modern Blender 4.4 APIs
    
    Args:
        star_count: Number of stars to generate
        distribution: Distribution pattern ('random', 'clustered', 'spiral')
        scale: Scale factor for star sizes
        **kwargs: Additional parameters
    
    Returns:
        dict: Created objects and metadata
    """
    # Use modern operator
    bpy.ops.astronomical.create_star_field(
        count=star_count,
        distribution=distribution,
        scale=scale
    )
    
    return {
        "objects": [obj for obj in bpy.context.scene.objects if obj.name.startswith("Star")],
        "type": "stellar_field",
        "metadata": {"count": star_count, "distribution": distribution, "scale": scale}
    }

def create_galaxy(galaxy_type="spiral", arms=4, bulge_size=1.0, **kwargs):
    """
    Create a galaxy using modern geometry nodes
    
    Args:
        galaxy_type: Type of galaxy ('spiral', 'elliptical', 'irregular')
        arms: Number of spiral arms (for spiral galaxies)
        bulge_size: Size of galactic bulge
        **kwargs: Additional parameters
    
    Returns:
        dict: Created objects and metadata
    """
    bpy.ops.astronomical.create_galaxy(
        galaxy_type=galaxy_type,
        arms=arms,
        bulge_size=bulge_size
    )
    
    return {
        "objects": [obj for obj in bpy.context.scene.objects if obj.name.startswith("Galaxy")],
        "type": "galaxy",
        "metadata": {"galaxy_type": galaxy_type, "arms": arms, "bulge_size": bulge_size}
    }

def setup_astronomical_scene(preset="scientific", **kwargs):
    """
    Setup complete astronomical scene with modern Blender 4.4 features
    
    Args:
        preset: Scene preset ('scientific', 'cinematic', 'publication')
        **kwargs: Override parameters
    
    Returns:
        dict: Scene configuration
    """
    # Setup scene with modern APIs
    setup_scene()
    setup_camera()
    setup_lighting()
    
    # Configure rendering engine (use EEVEE Next for Blender 4.4)
    bpy.context.scene.render.engine = 'BLENDER_EEVEE_NEXT'
    setup_rendering()
    
    return {
        "preset": preset,
        "render_engine": bpy.context.scene.render.engine,
        "scene": bpy.context.scene.name
    }

# ==========================================
# Module Exports
# ==========================================

__all__ = [
    # Registration
    "register",
    "unregister",
    # Main API
    "create_stellar_field", 
    "create_galaxy",
    "setup_astronomical_scene",
    # Core functions
    "setup_camera",
    "setup_lighting",
    "setup_rendering", 
    "setup_scene",
    "render_astronomical_scene",
    # Data utilities
    "STELLAR_CLASSIFICATION",
    "GALAXY_TYPES",
    "HR_DIAGRAM_PARAMS",
    "get_stellar_data",
    "get_galaxy_config",
    "create_sample_stellar_data",
    "validate_hr_diagram_data",
]

# ==========================================
# Entry Point
# ==========================================

if __name__ == "__main__":
    register()
