#Stereo Reconstruction

## Reconstruction:
The idea here is to capture the shape and appearance of objects from multiple images. The traditional stereo reconstruction pipeline consists of two main stages:
* Camera Calibration: Involves retrieval of the external camera pose and orientation as well as the internal parameters.
* Dense Matching: Establishing dense correspondences allows for a full 3D reconstruction and is generally done by global discrete or continuous optimization on the rectified images.

The estimation of the epipolar geometry can benefit from a large amount of accurate pixel correspondences, such as provided by optical flow methods.

## What does it Do:
Given two stereo image pairs, the aim is to reconstruct 3D models using simple linear triangulation methods. And finally to project the 3D points back into the images and calculate and report the mean square errors for each algorithm used and the two image pairs.

## Conclusion
* it is assumed that the cameras are either completely calibrated, or that the epipolar geometry of the image pair is known a priory. In such case, the image correspondence problem reduces to a one dimensional search along the known epipolar lines.
* These algorithms do not solve the problems completely, there still might be some occlusion and illumination changes.
* A best solution requires the definition and minimization of a suitable cost function and Optimal triangulation achieves that to some extent.
* In affine and projective reconstruction if there is no meaningful metric information about the object space. It is desirable to find a triangulation method that is invariant to projective transformations of space.
* An optimal solution requires avoiding a non-linear minimization of a cost function and hence the root of a sixth-degree polynomial.
* The homogeneous linear method described above often provides acceptable results but it is always recommended to use triangulation that is invariant to the projective frame of the cameras, and minimizes a geometric image error.
* The computational cost increases with an increase in the number of views.

### Usage: 
The folder contains all the inputs, just run the main.m in MATLAB, select the image pair that you wish to execute and all the required steps will be executed.
