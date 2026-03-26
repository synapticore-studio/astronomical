# Astronomical

[![ko-fi](https://ko-fi.com/img/githubbutton_sm.svg)](https://ko-fi.com/N4N71WOHZ3)

Professional astronomical visualization tools for Blender 4.4+

![Blender](https://img.shields.io/badge/Blender-4.4+-orange.svg)
![License](https://img.shields.io/badge/License-GPL--3.0-blue.svg)
![Version](https://img.shields.io/badge/version-2.0.0-green.svg)

[![GitHub](https://img.shields.io/badge/GitHub-Repository-black?logo=github)](https://github.com/bjoernbethge/astronomical)
[![Issues](https://img.shields.io/badge/GitHub-Issues-blue?logo=github)](https://github.com/bjoernbethge/astronomical/issues)

## Features

🌌 **Stellar Field Generation** - Create realistic star fields with proper stellar classification
🗺️ **Coordinate Systems** - ICRS, Galactic, Ecliptic transformations with accurate astronomical positioning
📊 **Catalog Integration** - Gaia, Hipparcos, Messier catalogs with real astronomical data
🎬 **Animation Tools** - Orbital mechanics, proper motion, and time-based astronomical phenomena
🎨 **Custom Node Groups** - Geometry nodes, shader nodes, and compositing nodes for astronomical rendering
💫 **Galaxy Morphology** - Realistic galaxy visualization with proper structure and classification
🌐 **Cosmic Web Simulation** - Large-scale structure simulation with filaments and voids
🔬 **Scientific Materials** - Physically-based materials for astronomical objects

## Installation

### Method 1: Blender Extensions Platform (Recommended)

1. Download the latest release from [GitHub Releases](https://github.com/bjoernbethge/astronomical/releases)
2. In Blender: **Edit → Preferences → Get Extensions**
3. Click the **Install from Disk** button (top-right)
4. Select the downloaded `.zip` file
5. Enable **Astronomical Blender Extension** in the extensions list

### Method 2: Manual Installation

1. Clone or download this repository
2. Run the build script:
   ```bash
   # Linux/macOS
   ./build.sh

   # Windows
   build.bat
   ```
3. Install the generated `.zip` file in Blender (see Method 1, steps 2-5)

## Quick Start

### Basic Usage

Once installed, access the Astronomical tools from the 3D Viewport sidebar:

1. Press `N` to open the sidebar
2. Navigate to the **Astronomical** tab
3. Choose from available operators:
   - **Create Stellar Field** - Generate realistic star fields
   - **Add Galaxy** - Create galaxy structures
   - **Setup Cosmic Web** - Generate large-scale structure
   - **Import Catalog** - Load real astronomical data

### Python API

```python
import bpy
from astronomical import (
    setup_scene,
    setup_camera,
    setup_lighting,
    render_astronomical_scene,
    get_stellar_data,
    create_sample_stellar_data
)

# Setup astronomical scene
setup_scene(bpy.context.scene)
setup_camera(bpy.context.scene)
setup_lighting(bpy.context.scene)

# Create stellar field with classification
stellar_data = create_sample_stellar_data(num_stars=1000)
bpy.ops.astronomical.create_stellar_field(count=1000)

# Render the scene
render_astronomical_scene(bpy.context.scene, output_path="/tmp/astronomy.png")
```

## 🎨 Visual Examples

### Stellar Field Visualization
Create realistic star fields with proper stellar classification and HR diagram positioning:

```python
# Create stellar field with spectral types
bpy.ops.astronomical.create_stellar_field(
    count=5000,
    distribution='realistic',  # Based on HR diagram
    magnitude_range=(0, 15)    # Apparent magnitude range
)
```

**Features:**
- Spectral classification (O, B, A, F, G, K, M types)
- HR diagram-based distribution
- Realistic color temperature and luminosity
- Proper motion animation support

### Galaxy Morphology
Visualize different galaxy types with scientifically accurate structures:

```python
# Create elliptical galaxy
bpy.ops.astronomical.add_galaxy(
    galaxy_type='elliptical',
    bulge_size=2.0,
    disk_size=0.0
)

# Create spiral galaxy with arms
bpy.ops.astronomical.add_galaxy(
    galaxy_type='spiral',
    num_arms=2,
    arm_width=0.5,
    rotation=45.0
)
```

**Galaxy Types:**
- Elliptical (E0-E7)
- Spiral (Sa, Sb, Sc)
- Barred Spiral (SBa, SBb, SBc)
- Irregular (Irr)

### Cosmic Web Simulation
Generate large-scale structure with filaments and voids:

```python
# Setup cosmic web simulation
bpy.ops.astronomical.setup_cosmic_web(
    grid_size=100,
    filament_density=0.3,
    void_probability=0.2
)
```

**Features:**
- Filamentary structure generation
- Void region simulation
- Galaxy cluster positioning
- Matter density visualization

### Coordinate System Transformations
Work with multiple astronomical coordinate systems:

```python
from astronomical.utilities.astronomical_data import (
    transform_coordinates,
    COORDINATE_SYSTEMS
)

# Transform from ICRS to Galactic coordinates
galactic_coords = transform_coordinates(
    ra=123.45, dec=67.89,
    from_system='icrs',
    to_system='galactic'
)
```

**Supported Systems:**
- ICRS (International Celestial Reference System)
- Galactic (Galactic coordinates)
- Ecliptic (Ecliptic coordinates)

> 📸 **Screenshots coming soon** - Beautiful astronomical visualizations!

> **Note**: High-quality screenshots and demo renders will be added soon. Track progress in [GitHub Issues](https://github.com/bjoernbethge/astronomical/issues).

## Requirements

- **Blender**: 4.4.0 or higher
- **Python**: 3.11+ (bundled with Blender)
- **Platforms**: Windows, macOS, Linux (x64 and ARM64)

### Optional Dependencies

For advanced features, the extension can use:
- **NumPy**: Numerical computations (bundled with Blender)
- **AstroPy**: Astronomical calculations (vendored)
- **SciPy**: Scientific algorithms (optional)

## Blender Compatibility

| Blender Version | Status | Notes |
|----------------|--------|-------|
| 4.4.x          | ✅ Full Support | Recommended version |
| 4.5.x          | ✅ Full Support | Latest features |
| 4.3.x          | ⚠️ Limited | Some features unavailable |
| 4.2.x or older | ❌ Not Supported | Use legacy addon version |

## Node Groups

The extension includes custom node groups for advanced astronomical rendering:

### Geometry Nodes
- **Stellar Point Cloud** - Generate point-based star fields
- **Galaxy Spiral Arms** - Procedural spiral galaxy generation
- **Cosmic Web Filaments** - Large-scale structure filaments

### Shader Nodes
- **Star Surface Shader** - Realistic stellar surface materials
- **Nebula Volume** - Volumetric nebula rendering
- **Galaxy Disk Material** - Galaxy disk with dust lanes

### Compositing Nodes
- **Star Diffraction** - Diffraction spike generation
- **Atmospheric Scattering** - Atmospheric effects for Earth views
- **Light Pollution** - Urban light pollution simulation

## Development

### Building from Source

```bash
# Clone repository
git clone https://github.com/bjoernbethge/astronomical.git
cd astronomical

# Build extension package
./build.sh  # Linux/macOS
# or
build.bat   # Windows

# The output will be in the root directory: astronomical_blender_extension-2.0.0.zip
```

### Project Structure

```
astronomical/
├── __init__.py              # Main extension entry point
├── blender_manifest.toml    # Extension metadata (Blender 4.4+)
├── core/                    # Core functionality
│   ├── camera.py           # Camera setup and positioning
│   ├── lighting.py         # Lighting systems
│   ├── rendering.py        # Render configuration
│   └── scene.py            # Scene setup utilities
├── nodes/                   # Custom node groups
│   ├── geometry/           # Geometry nodes
│   ├── shader/             # Shader nodes
│   └── compositing/        # Compositing nodes
├── operators/              # Blender operators
│   ├── stellar_field.py   # Stellar field generation
│   ├── galaxy.py          # Galaxy creation
│   └── cosmic_web.py      # Cosmic web simulation
└── utilities/              # Utility modules
    ├── astronomical_data.py  # Astronomical data and constants
    └── transforms.py         # Coordinate transformations
```

### Testing

The extension includes a test suite:

```bash
# Run tests (requires Blender Python)
blender --background --python-console
>>> import astronomical
>>> astronomical.test()
```

## Contributing

Contributions are welcome! Please follow these guidelines:

1. **Code Style**: Follow PEP 8 and Blender addon conventions
2. **Modularity**: Keep functions small and focused (DRY principle)
3. **Documentation**: Document all public functions and classes
4. **Testing**: Add tests for new features
5. **Pull Requests**: Create feature branches and submit PRs

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

## Use Cases

### Scientific Visualization
- Research paper illustrations
- Educational materials
- Planetarium content
- Documentary visualizations

### Artistic Projects
- Science fiction environments
- Space-themed artworks
- Procedural starscapes
- Cosmic landscapes

### Data Visualization
- Astronomical survey data
- N-body simulation results
- Observational data from telescopes
- Statistical distributions in 3D

## Roadmap

### Upcoming Features
- [ ] Real-time catalog data import from online databases
- [ ] Advanced animation tools for orbital mechanics
- [ ] Integration with external astronomical software
- [ ] VR support for immersive space exploration
- [ ] Machine learning-based galaxy generation
- [ ] Enhanced nebula rendering with spectral data

Track feature requests and vote on priorities in [GitHub Issues](https://github.com/bjoernbethge/astronomical/issues).

## License

**GPL-3.0-or-later** - see [LICENSE](LICENSE) for details.

This extension is free and open-source software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.

## Credits

Built with:
- [Blender](https://www.blender.org/) - 3D creation suite (4.4+)
- [AstroPy](https://www.astropy.org/) - Astronomical calculations (vendored)
- [NumPy](https://numpy.org/) - Numerical computations (bundled)

Special thanks to the astronomical visualization community and contributors.

## Related Links

**Project:**
- [GitHub Repository](https://github.com/bjoernbethge/astronomical) - Source code and issues
- [Report Issues](https://github.com/bjoernbethge/astronomical/issues) - Bug reports and feature requests
- [Releases](https://github.com/bjoernbethge/astronomical/releases) - Download latest version

**Documentation & Resources:**
- [Blender Documentation](https://docs.blender.org/) - Official Blender docs
- [Blender 4.4 Release Notes](https://www.blender.org/download/releases/4-4/) - What's new in Blender 4.4
- [Blender Python API](https://docs.blender.org/api/current/) - bpy API reference
- [AstroPy Documentation](https://docs.astropy.org/) - Astronomical calculations

**Astronomical Data Sources:**
- [Gaia Archive](https://gea.esac.esa.int/archive/) - European Space Agency star catalog
- [SDSS](https://www.sdss.org/) - Sloan Digital Sky Survey
- [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/) - Exoplanet data
- [Messier Catalog](https://en.wikipedia.org/wiki/Messier_object) - Deep sky objects

---

**Astronomical** - Professional astronomical visualization for Blender 🌌✨
