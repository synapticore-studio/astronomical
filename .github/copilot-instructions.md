# Astronomical Blender Extension - Copilot Instructions

## Project Overview

This is a professional astronomical visualization extension for Blender 4.4+, providing tools for creating realistic stellar fields, galaxy morphology, cosmic web simulations, and scientific astronomical rendering.

## Tech Stack

- **Blender**: 4.4.0 or higher (target platform)
- **Python**: 3.11+ (bundled with Blender)
- **Core Dependencies**: NumPy (bundled), AstroPy (vendored), SciPy (optional)
- **Architecture**: Blender Extension (modern manifest-based, not legacy addon)
- **Platforms**: Windows, macOS, Linux (x64 and ARM64)

## Project Structure

```
astronomical/
├── __init__.py              # Main extension entry point and registration
├── blender_manifest.toml    # Extension metadata (Blender 4.4+ format)
├── core/                    # Core functionality modules
│   ├── camera.py           # Camera setup and positioning
│   ├── lighting.py         # Lighting systems
│   ├── rendering.py        # Render configuration
│   ├── scene.py            # Scene setup utilities
│   └── data.py             # Data management
├── nodes/                   # Custom node groups
│   ├── geometry/           # Geometry nodes for procedural generation
│   ├── shader/             # Shader nodes for materials
│   └── compositing/        # Compositing nodes for post-processing
├── operators/              # Blender operators (user-facing operations)
│   ├── star.py            # Stellar object creation
│   ├── galaxy.py          # Galaxy generation
│   ├── pointcloud.py      # Point cloud operations
│   └── ...                # Other astronomical operators
└── utilities/              # Utility modules
    ├── astronomical_data.py   # Astronomical constants and data
    ├── astropy_integration.py # AstroPy wrapper functions
    └── ...                    # Other utilities
```

## Code Style and Conventions

### Python Style
- **Follow PEP 8** strictly for all Python code
- **Follow Blender addon conventions** for operator and property naming
- Use **type hints** for all function parameters and returns
- Use **docstrings** for all public classes, functions, and modules (Google style)
- Keep functions **small and focused** (DRY principle)
- Use **meaningful variable names** (e.g., `spectral_class`, not `sc`)

### Blender-Specific Conventions
- **Operator naming**: Use `ALBPY_OT_` prefix for operator classes (e.g., `ALBPY_OT_CreateStar`)
- **Property naming**: Use descriptive names with proper units in descriptions
- **bl_idname format**: Use `albpy.operation_name` (lowercase, underscores)
- **Registration**: All operators, panels, and properties must be properly registered/unregistered
- **Modern APIs**: Use Blender 4.4+ APIs (e.g., `bpy.utils.expose_bundled_modules()`)
- **Node Interface API**: Use modern node interface API, not legacy socket creation

### Module Organization
- **All imports at top of file** (standard library, third-party, local)
- **Group related functionality** in logical modules (camera, lighting, rendering)
- **Avoid circular imports** - use proper dependency hierarchy
- **Vendor dependencies** in a `vendor/` directory when needed

### Error Handling
- Use **try-except blocks** for operations that may fail
- **Log errors properly** using `print()` statements (Blender console)
- **Provide user feedback** through `self.report()` in operators
- **Validate inputs** before processing

## Build and Testing

### Building
```bash
# Linux/macOS
./build.sh

# Windows
build.bat

# Quick build (development)
python3 quick_build.py

# Full build
python3 build_extension.py

# Clean build
python3 build_extension.py --clean
```

Output: `astronomical_blender_extension-2.0.0.zip`

### Testing
```bash
# Run tests in Blender Python console
blender --background --python-console
>>> import astronomical
>>> astronomical.test()
```

### Installation Testing
1. Build the extension
2. Install in Blender: Edit → Preferences → Get Extensions → Install from Disk
3. Enable the extension in the extensions list
4. Test in 3D Viewport sidebar (press N, find "Astronomical" tab)

## Development Workflow

### Adding New Features
1. **Identify the correct module** (core, operators, nodes, utilities)
2. **Follow existing patterns** in similar files
3. **Create operators** for user-facing operations
4. **Add properties** using `bpy.props` with proper descriptions
5. **Register/unregister** all new classes properly
6. **Document** with docstrings and comments where needed
7. **Test** in Blender before committing

### Modifying Existing Code
1. **Understand the existing implementation** before making changes
2. **Maintain backward compatibility** where possible
3. **Update docstrings** if behavior changes
4. **Test thoroughly** to ensure no regressions
5. **Keep changes minimal** and focused

### Coordinate Systems
When working with astronomical coordinates:
- **ICRS**: International Celestial Reference System (standard)
- **Galactic**: Galactic coordinates (l, b)
- **Ecliptic**: Ecliptic coordinates
- Use `utilities/astronomical_data.py` for transformations

### Stellar Classification
When working with stars:
- Use **spectral classes**: O, B, A, F, G, K, M (hottest to coolest)
- Include **temperature** in Kelvin (e.g., Sun = 5778K)
- Include **luminosity** in solar luminosities
- Use **HR diagram parameters** from `utilities/astronomical_data.py`

## Boundaries and Restrictions

### Do NOT
- ❌ **Do not modify** `blender_manifest.toml` without understanding Blender extension format
- ❌ **Do not use** deprecated Blender APIs (pre-4.4 patterns)
- ❌ **Do not break** existing operator functionality
- ❌ **Do not commit** temporary files, `__pycache__`, `.pyc` files
- ❌ **Do not use** `any` type hints - be specific with types
- ❌ **Do not add** external dependencies without vendoring them properly
- ❌ **Do not modify** core registration logic without careful testing
- ❌ **Do not remove** existing tests or functionality without justification

### Always DO
- ✅ **Do follow** PEP 8 and Blender addon conventions
- ✅ **Do use** type hints for all functions
- ✅ **Do document** public APIs with docstrings
- ✅ **Do test** in Blender before committing
- ✅ **Do validate** user inputs in operators
- ✅ **Do use** proper error handling
- ✅ **Do keep** functions focused and modular
- ✅ **Do check** Blender version compatibility (4.4+)

### Dependencies
- **NumPy**: Already bundled with Blender - use freely
- **AstroPy**: Vendored in project - use via vendor path
- **New dependencies**: Must be vendored or confirmed available in Blender Python
- **Check availability** before using any external library

## Scientific Accuracy

This extension aims for scientific accuracy in astronomical visualization:
- **Use real astronomical constants** from `utilities/astronomical_data.py`
- **Implement proper physics** (inverse square law, proper motion, etc.)
- **Follow astronomical conventions** (coordinate systems, units, classifications)
- **Cite sources** when implementing algorithms based on papers/datasets
- **Validate data ranges** (e.g., temperature ranges for spectral classes)

## UI/UX Guidelines

### Operator Design
- **Clear descriptions**: Every operator needs a clear `bl_description`
- **Logical property grouping**: Group related properties together
- **Sensible defaults**: Use scientifically reasonable default values
- **Unit labels**: Always include units in property descriptions (K, solar radii, etc.)
- **Undo support**: Use `'UNDO'` in `bl_options` where appropriate

### Panel Layout
- **Logical organization**: Group related tools in the same panel
- **Clear labels**: Use descriptive labels for all UI elements
- **Tooltips**: Provide helpful descriptions in property definitions
- **Progressive disclosure**: Show advanced options only when needed

## Common Patterns

### Creating Objects
```python
def execute(self, context):
    # Create mesh data
    mesh = bpy.data.meshes.new(name="ObjectName")
    obj = bpy.data.objects.new(name="ObjectName", object_data=mesh)
    
    # Link to scene
    context.collection.objects.link(obj)
    
    # Set as active
    context.view_layer.objects.active = obj
    obj.select_set(True)
    
    return {'FINISHED'}
```

### Using Properties
```python
my_float: FloatProperty(
    name="Property Name",
    description="Clear description with units (if applicable)",
    default=1.0,
    min=0.0,
    max=10.0,
    step=0.1,
    precision=2
)
```

### Error Handling in Operators
```python
def execute(self, context):
    try:
        # Main operation
        result = perform_operation()
        self.report({'INFO'}, "Operation completed successfully")
        return {'FINISHED'}
    except Exception as e:
        self.report({'ERROR'}, f"Operation failed: {str(e)}")
        return {'CANCELLED'}
```

## Resources

- **Blender Python API**: https://docs.blender.org/api/current/
- **AstroPy Documentation**: https://docs.astropy.org/
- **Gaia Archive**: https://gea.esac.esa.int/archive/
- **Project Repository**: https://github.com/bjoernbethge/astronomical

## Questions?

If you're unsure about any convention or pattern:
1. Check existing similar code in the repository
2. Refer to Blender documentation for API usage
3. Follow the patterns established in this file
4. When in doubt, prefer clarity and simplicity over cleverness
