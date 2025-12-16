
````markdown
# Tracking Coordinate Transforms


Coordinate system and reference-frame transformations for tracking systems.

This repository contains an initial implementation of the coordinate
transformation layer on the Earth ellipsoid. Input data may be either in a coordinate system (azimuth, range, height) or geodetic coordinates (latitude, longitude, height).

The module transforms these representations into several local 2D Cartesian
coordinate systems, suitable for planar filtering. 

A detailed mathematical description of the underlying coordinate system
will be added later.


## Features
- Transform radar coordinates (azimuth, range, height) into geodetic coordinates (latitude, longitude, height) on the Earth Ellipsoid

Transform geodetic coordinates on the Earth Ellipsoid (latitude, longitude, height) into a local 2D Cartesian system (x,y) using three methods:
- Custom internal coordinate system based on length-preserving meridians and parallels
- Stereo projection
- Scaled geodetic system
- ENU system

Domain-agnostic: Applicable to air traffic, satellites, ships, ground vehicles, or any moving objects.


## Installation
Clone the repository:
```bash
git clone https://github.com/maksimenko421/tracking-coordinate-transforms.git
cd tracking-coordinate-transforms
````

## Configuration

Edit `general_settings.py` to configure Earth ellipsoid


## Usage

Example Python usage:

run files starting with test_*

The graph_on = True / False parameter enables/disables the round-trip error graphs.


## Mathematical background

A detailed mathematical description will be added in a future update.


## License

MIT License

Copyright (c) 2025 Irina Maksimenko

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE
# tracking-coordinate-transforms
