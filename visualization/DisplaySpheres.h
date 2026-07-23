#ifndef DISPLAYSPHERES_INCLUDED
#define DISPLAYSPHERES_INCLUDED

#include <vector>
#include <functional>
//editted by Jaeuk
#include "PeriodicCellList.h"

struct sphere
{
	double x, y, z;//coordinates, on the order of 1
	               //the rotation center is (0, 0, 0) 
	double radius;

	double red, green, blue;//color, from 0 to 1, 0:dark, 1:bright
	double transparency;//from 0 to 1, 0:invisible, 1:opaque
};
struct line
{
	double x1, y1, z1;
	double x2, y2, z2;

	double width;// on the order of 1

	double red, green, blue;//color, from 0 to 1, 0:dark, 1:bright
	double transparency;//from 0 to 1, 0:invisible, 1:opaque
};

//should only be called inside the Extra2DStuff parameters of DisplaySpheres functions.
/** \brief Draw 2D overlay text; only valid inside the Extra2DStuff callback of DisplaySpheres.
 *  \param text  Text string to draw.
 *  \param x     X position of the text in window coordinates.
 *  \param y     Y position of the text in window coordinates. */
void DisplayText(const std::string & text, double x=0.0, double y=0.0);

// Display all spheres in the vector
// Use mouse to drag spheres
// Use keyboard arrows to move the camera
// Use keys '+' and '-' to increase or decrease the precision. Spheres looks bad if precision is too low. Computer is laggy if precision is too high.
// Use key 'q' to return from this function
// Use key 'o' to switch to orthogonal projection
// Use key 'p' to switch to perspective projection (default)
// Use key 'v' to write a Pov-ray source file of current configuration and view angle
// Don't click the Close Window icon in the upper right corner of the title bar. That will kill the program.
// In order for this function to work, define USE_OPENGL when compiling, and link with OpenGL, GLU, and GLUT
/** \brief Display all spheres and lines in an interactive OpenGL window.
 *  \param VSpheres      Spheres to render.
 *  \param VLines        Lines to render.
 *  \param Extra2DStuff  Callback receiving (window width, height) to draw extra 2D overlay content. */
void DisplaySpheres(const std::vector<sphere> & VSpheres, const std::vector<line> & VLines, std::function<void(size_t, size_t)> Extra2DStuff=[](size_t, size_t)->void{});


//Similar to DisplaySpheres, except renders a movie
// Use key 'a' to accelerate the movie
// Use key 's' to deccelerate the movie
// Use key 'd' to toggle pause/play (pause at the beginning)
// Use key 'r' to toggle forward/backward (forward by default)
/** \brief Display an animated movie of sphere/line frames in an interactive OpenGL window.
 *  \param GetFrameFunc  Callback that fills the spheres and lines for a given frame number.
 *  \param NumFrames     Total number of frames in the movie. */
void DisplaySpheres_Movie(std::function<void(std::vector<sphere> & VSpheres, std::vector<line> & VLines, long FrameNumber)> GetFrameFunc, long NumFrames);

//transferred from PeriodicCellList.h
/** \brief Display a particle configuration as 3D spheres (with optional bonds) in an OpenGL window.
 *  \param a           Configuration to display (1D, 2D, or 3D).
 *  \param BondLength  Maximum neighbor distance to draw bonds; 0 draws no bonds.
 *  \param Radius      Sphere radius; 0 selects an automatic radius. */
void DisplayConfiguration(const Configuration & a, double BondLength = 0.0, double Radius = 0.0);

/** \brief Display a configuration pack as an animated 3D movie in an OpenGL window.
 *  \param pk          Configuration pack providing one configuration per frame.
 *  \param BondLength  Maximum neighbor distance to draw bonds; 0 draws no bonds.
 *  \param Radius      Sphere radius; 0 selects an automatic radius. */
void Display3DConfigurationMovie(ConfigurationPack & pk, double BondLength = 0.0, double Radius = 0.0);

/** \brief Display the 3D structure factor of a configuration as colour-coded spheres.
 *  \param a              3D configuration to analyze.
 *  \param KCutOff        Maximum wavevector magnitude to include.
 *  \param KSphereRadius  Radius drawn for each k-point sphere; 0 selects an automatic radius. */
void Display3DStructureFactor(const Configuration& a, double KCutOff, double KSphereRadius = 0.0);

#endif