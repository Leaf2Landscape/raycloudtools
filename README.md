## Ray Cloud Tools
A set of command line tools for processing ray clouds, together with an associated C++ library.

If you publish articles based on this tool set, please cite
```
Lowe, Thomas, et al. "Canopy Density Estimation in Perennial Horticulture Crops Using 3D Spinning LiDAR SLAM." arXiv preprint arXiv:2007.15652 (2020).
```
### Paper (bibtex)
```
@article{lowe2020canopy,
  title={Canopy Density Estimation in Perennial Horticulture Crops Using 3D Spinning LiDAR SLAM},
  author={Lowe, Thomas and Moghadam, Peyman and Edwards, Everard and Williams, Jason},
  journal={arXiv preprint arXiv:2007.15652},
  year={2020}
}
```
 

## Build:
```console
mkdir build
cd build
cmake ..
```

To access the tools from anywhere, place in your ~/bashrc:
```console
  export PATH=$PATH:'source code path'/raycloudtools/build/bin
```

*Dependencies:*

Eigen: sudo apt-get install libeigen3-dev

LibNabo: 
```console
git clone https://github.com/ethz-asl/libnabo.git
git checkout tags/1.0.7
``` 
then follow build and install instructions in its README.md.

## Examples:

**rayimport forest.laz forest_traj.ply** &nbsp;&nbsp;&nbsp; Import point cloud and trajectory to a single raycloud file: forest.ply

**raycreate room 1** &nbsp;&nbsp;&nbsp; Generate a single room with a window and door, using random seed 1.
<p align="center">
<img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room1.png?at=refs%2Fheads%2Fmaster"/>
  <img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room3.png?at=refs%2Fheads%2Fmaster"/>
  <img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room2.png?at=refs%2Fheads%2Fmaster"/>
</p>
&nbsp;&nbsp;&nbsp; You can visualise the rays in meshlab with Render | Show Vertex Normals. The ray lengths need to be scaled: Tools | Options | NormalLength roughly 0.025 (smaller for larger clouds)

**raydecimate room.ply 10 cm** &nbsp;&nbsp;&nbsp; Spatially decimate cloud to one point every cubic 10 cm.

<p align="center"><img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room_decimated.png?at=refs%2Fheads%2Fmaster"/></p>

**raytranslate room.ply 3 0 0** &nbsp;&nbsp;&nbsp; Translate the ray cloud 3 metres along the x axis.

<p align="center"><img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room_translate.png?at=refs%2Fheads%2Fmaster"/></p>

**rayrotate room.ply 0 0 30** &nbsp;&nbsp;&nbsp; Rotate the ray cloud 30 degrees around the z axis.

<p align="center"><img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room_rotate.png?at=refs%2Fheads%2Fmaster"/></p>

**raydenoise room.ply 10 cm** &nbsp;&nbsp;&nbsp; Remove rays with isolated end points more than 10 cm from any other, not including unbounded rays.

<p align="center">
<img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room_denoise1.png?at=refs%2Fheads%2Fmaster"/>
<img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room_denoise2.png?at=refs%2Fheads%2Fmaster"/>

**raysmooth room.ply** &nbsp;&nbsp;&nbsp; Move ray end points onto the nearest surface, to smooth the resulting cloud.

<p align="center">
<img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room_smooth1.png?at=refs%2Fheads%2Fmaster"/>
<img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room_smooth2.png?at=refs%2Fheads%2Fmaster"/>
</p>

**raytransients min room.ply 2 rays** &nbsp;&nbsp;&nbsp; Segment out moving or moved objects during the scan, when matter has been re-observed as missing by 2 or more rays. 

&nbsp;&nbsp;&nbsp; Leaving the ***minimum*** of geometry when transient.

&nbsp;&nbsp;&nbsp; In this raycloud the table and cupboard appear only after the empty room has been scanned for several seconds, so we can isolate these transient objects.

<p align="center">
<img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room_transients1.png?at=refs%2Fheads%2Fmaster"/>
<img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room_transients2.png?at=refs%2Fheads%2Fmaster"/>
<img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room_transients3.png?at=refs%2Fheads%2Fmaster"/>
</p>

&nbsp;&nbsp;&nbsp; Left: original cloud. Middle: the fixed (untransient) raycloud. Right: the remaining transient rays are also saved.

**raycombine all room.ply room2.ply** &nbsp;&nbsp;&nbsp; Combine room and its transformed version together, keeping ***all*** rays.

<p align="center"><img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room_combined_all.png?at=refs%2Fheads%2Fmaster"/></p>

**raycombine min room.ply room2.ply 1 rays** &nbsp;&nbsp;&nbsp; Combine the two ray clouds keeping only the ***minimum*** of geometry where there is a difference. 

&nbsp;&nbsp;&nbsp; This is a form of union of the two volumes. 

<p align="center"><img img width="320" src="https://raw.githubusercontent.com/csiro-robotics/raycloudtools/main/pics/room_combined_min.png?at=refs%2Fheads%2Fmaster"/></p>

**rayalign room.ply room2.ply** &nbsp;&nbsp;&nbsp; Aligns room onto room2, allowing for a small about of non-rigidity 

*Optional build dependencies:*

For rayconvert to work from .laz files:
* git clone https://github.com/LASzip/LASzip.git, then git checkout tags/2.0.1, then mkdir build, cd build, cmake .., make, sudo make install. 
* git clone https://github.com/libLAS/libLAS.git, then mkdir build, cd build, cmake .. -DWITH_LAS_ZIP=ON, make, sudo make install
* in raycloudtools/build: cmake .. -DWITH_LAS=ON  (or ccmake .. to turn on/off WITH_LAS)

For raywrap:

* git clone http://github.com/qhull/qhull.git, git checkout tags/v7.3.2
* In qhull: mkdir build, cd build, cmake .. -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true, make, sudo make install. 
* in raycloudtools/build: cmake .. -DWITH_QHULL=ON (or ccmake .. to turn on/off WITH_QHULL)

## Unit Tests

Unit tests must be enabled at build time before running. To build with unit tests, the CMake variable `RAYCLOUD_BUILD_TESTS` must be `ON`. This can be done in the initial project configuration by running the following command from the `build` directory: `cmake  -DRAYCLOUD_BUILD_TESTS=ON ..`

Unit tests may then be run directly or using `CTest`.

### Running using CTest

To run using CTest:

* Change into the `build` directory
* Run `ctest .`

On some platforms it may be necessary to specify the build configuration to test. For example, the `Release` build may be tested using the modified command `ctest . -C Release`.

### Directly invoking the tests

When directly invoking the unit tests, is important that the tests are run from the directory to which the raycloud tools executables are built. To invoke the tests directly:

* Change into the `build` directory
* Change into the `bin/` directory
* Run `./raytest`

## Acknowledgements
This research was supported by funding from CSIRO's Data61, Land and Water, Wine Australia, and the Department of Agriculture's Rural R&D for Profit program. The authors gratefully acknowledge the support of these groups, which has helped in making this library possible. 


## Notes

### Dependencies

> To be converted into instructions.

* liblas apt packages `liblas-dev liblas-c-dev`
* qhull: apt packages do not include cpp libraries (Ubuntu 18.04)
* ros: requires `source /opt/ros/<distro>/setup.bash` before configuring cmake

### Algorithmic

* Need to support different reference frames in debug visualisation.
  * Check current assumptions: assuming Z up, right handed?
