#!/usr/bin/env python3
"""
Astronomical Blender Extension Build Script
==========================================

Automated build system for the AstroPhot Blender Extension.
Handles dependency packaging, wheel bundling, and extension creation.

Usage:
    python build_extension.py [--clean] [--dev] [--wheels-only]
    
Options:
    --clean      Clean build directories before building
    --dev        Development build (skip dependency packaging)
    --wheels-only Only download and package