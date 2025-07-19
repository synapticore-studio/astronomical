#!/bin/bash
# Simple build script for macOS/Linux
# Following project standards - direct approach only

echo "Build Astronomical Extension"

# Check Python
if ! command -v python3 &> /dev/null; then
    echo "Python 3 not found"
    exit 1
fi

echo "Python version:"
python3 --version
echo

echo "Choose build:"
echo "1. Quick build"
echo "2. Full build" 
echo "3. Clean build"
echo

read -p "Choice (1-3): " choice

case $choice in
    1)
        python3 quick_build.py
        ;;
    2)
        python3 build_extension.py
        ;;
    3)
        python3 build_extension.py --clean
        ;;
    *)
        echo "Invalid choice"
        exit 1
        ;;
esac

echo
echo "Build complete"
