"""
Core Astronomical Extension Components
=====================================

Core functionality for scene setup, rendering, camera, and lighting.
Self-contained implementation without external dependencies.
"""

from .camera import setup_camera, create_orbital_camera, get_camera_info, calculate_fov
from .data import (
    get_coordinates_and_features, 
    list_available_surveys, 
    load_survey_data,
    get_data_statistics,
    validate_survey_data,
    create_sample_dataset
)
from .lighting import (
    setup_lighting, 
    create_star_lighting, 
    create_nebula_lighting,
    get_lighting_info,
    optimize_lighting_for_engine
)
from .rendering import (
    render_astronomical_scene, 
    setup_rendering,
    get_render_info,
    optimize_for_headless,
    setup_astronomical_post_processing
)
from .scene import (
    setup_scene,
    setup_astronomical_view_settings,
    reset_scene_to_defaults,
    get_scene_info,
    validate_scene_for_astronomy
)

__all__ = [
    # Camera functions
    "setup_camera",
    "create_orbital_camera", 
    "get_camera_info",
    "calculate_fov",
    # Data functions
    "get_coordinates_and_features",
    "list_available_surveys",
    "load_survey_data", 
    "get_data_statistics",
    "validate_survey_data",
    "create_sample_dataset",
    # Lighting functions
    "setup_lighting",
    "create_star_lighting",
    "create_nebula_lighting",
    "get_lighting_info",
    "optimize_lighting_for_engine",
    # Rendering functions
    "render_astronomical_scene",
    "setup_rendering",
    "get_render_info", 
    "optimize_for_headless",
    "setup_astronomical_post_processing",
    # Scene functions
    "setup_scene",
    "setup_astronomical_view_settings",
    "reset_scene_to_defaults",
    "get_scene_info",
    "validate_scene_for_astronomy",
]
