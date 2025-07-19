#!/usr/bin/env python3
"""
Quick Build - KISS Implementation
================================

Simple, fast build for development.
No complex features - just copy files and zip.
"""

import shutil
import zipfile
from pathlib import Path
import toml


def quick_build():
    """Simple development build"""
    print("Quick build starting...")
    
    project_root = Path(__file__).parent
    
    # Read manifest
    manifest = toml.load(project_root / "blender_manifest.toml")
    version = manifest["version"]
    extension_name = manifest["id"]
    
    # Simple build directory
    build_dir = project_root / "quick_build"
    if build_dir.exists():
        shutil.rmtree(build_dir)
    
    ext_dir = build_dir / extension_name
    ext_dir.mkdir(parents=True)
    
    # Copy core files only
    core_files = [
        "__init__.py",
        "blender_manifest.toml",
        "core/",
        "nodes/", 
        "operators/",
        "utilities/",
    ]
    
    for item in core_files:
        source = project_root / item
        dest = ext_dir / item
        
        if source.is_file():
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, dest)
        elif source.is_dir():
            shutil.copytree(source, dest, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    
    # Create simple ZIP
    zip_name = f"{extension_name}-{version}-dev.zip"
    zip_path = project_root / zip_name
    
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file_path in ext_dir.rglob("*"):
            if file_path.is_file():
                rel_path = file_path.relative_to(build_dir)
                zipf.write(file_path, rel_path)
    
    shutil.rmtree(build_dir)
    
    print(f"Quick build complete: {zip_path}")
    return zip_path


if __name__ == "__main__":
    quick_build()
