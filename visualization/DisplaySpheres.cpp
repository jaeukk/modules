// DisplaySpheres.cpp
// adapted from SphereWorld.cpp found in OpenGL SuperBible, by Richard S. Wright Jr.

#include "DisplaySpheres.h"

#ifdef USE_OPENGL
#include <etc.h>
#include "math3d.h"
#include "glframe.h"
#include <StructureFactor.h>	//For Display3DStructureFactor(...) function

//#include <GL/glut.h>
#include <GL/freeglut.h>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

static std::vector<GLFrame> spheres;
static GLFrame    frameCamera;
static GLFrame		AllSpheres;

static std::vector<sphere> Spheres;
static std::vector<size_t> RenderOrder;
static std::vector<line> Lines;
static std::function<void(std::vector<sphere> & VSpheres, std::vector<line> & VLines, long FrameNumber)> * pGetFrameFunc;
static std::function<void(size_t, size_t)> * pExtra2DStuffFunc;
// Light and material Data
static GLfloat fLightPos[4]   = { -100.0f, 100.0f, 100.0f, 1.0f };  // Point source
static GLfloat fNoLight[] = { 0.0f, 0.0f, 0.0f, 0.0f };
static GLfloat fAmbientLight[] = { 0.19225f, 0.19225f, 0.19225f, 1.0f };
static GLfloat fDiffuseLight[] = { 0.50754f, 0.50754f, 0.50754f, 1.0f };
static GLfloat fSpecularLight[] = { 0.508273f, 0.508273f, 0.508273f, 1.0f };
static GLfloat ClearColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };
static GLfloat shine=0.4;

static size_t Precision = 16;
static GLfloat fAspect;
static size_t WindowWidth=800;
static size_t WindowHeight=600;
static bool OrthoProjection=false;

static long FrameRate = 1000000; //0.1 second in Windows
static bool Playing = false;
static bool Reverse = false;
static long MaxFrame;
static long CurrentFrame;
static long LastFrameTime;

static int window_handle;

/** \brief Fetch frame data and rebuild the per-sphere GLFrames for the given frame. (CLAUDE)
 *  \param FrameNumber  Index of the frame to prepare. */
void PrepareFrame (long FrameNumber)
{
	(*pGetFrameFunc)(::Spheres, ::Lines, FrameNumber);
	CurrentFrame=FrameNumber;
	LastFrameTime=GetPreciseClock();

	spheres.resize(Spheres.size());
	for(int iSphere = 0; iSphere < Spheres.size(); iSphere++)
	{
		spheres[iSphere].SetOrigin(static_cast<float>(Spheres.at(iSphere).x), static_cast<float>(Spheres.at(iSphere).y), static_cast<float>(Spheres.at(iSphere).z));
	}
	RenderOrder.clear();
	for(size_t i=0; i<Spheres.size(); i++)
		RenderOrder.push_back(i);
}

//////////////////////////////////////////////////////////////////
// This function does any needed initialization on the rendering
// context.
/** \brief Perform any needed initialization on the rendering context. */
void SetupRC()
{
	// Grayish background
	glClearColor(ClearColor[0], ClearColor[1], ClearColor[2], ClearColor[3]);

	// Cull backs of polygons
	//glCullFace(GL_BACK);
	//glFrontFace(GL_CCW);
	//glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);

	// Setup light parameters
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, fNoLight);
	glLightfv(GL_LIGHT0, GL_AMBIENT, fAmbientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, fDiffuseLight);
	glLightfv(GL_LIGHT0, GL_SPECULAR, fSpecularLight);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	//glEnable(GL_POLYGON_SMOOTH);

	// Mostly use material tracking
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv(GL_FRONT, GL_SPECULAR, fSpecularLight);
	glMateriali(GL_FRONT, GL_SHININESS, shine*128);

	// place the spheres
	AllSpheres.SetOrigin(0.0, 0.0, -0.0);

	frameCamera.SetOrigin(0.0, 0.0, 3.0);

	PrepareFrame(0);
}

void DisplayText(const std::string & text, double x, double y)
{
	glColor3f(0.0, 0.0, 0.0);
	glRasterPos2f(x, y);
	for ( int i = 0; i < text.size(); ++i ) 
	{
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, text[i]);
	}
}

//save the figure as a povray script
/** \brief Save the current figure and view as a Pov-Ray source file. */
void SavePovRay(void)
{
	auto t = std::time(nullptr);

	struct tm * now = localtime( & t );
	std::stringstream ss;
	ss	<< (now->tm_year + 1900) 
		<< (now->tm_mon + 1) 
		<<  now->tm_mday 
		<< now->tm_hour 
		<< now->tm_min 
		<< now->tm_sec 
		<< ".pov";
	std::fstream ofile(ss.str(), std::fstream::out);
	M3DVector3f co, ca;
	ofile<<"#version 3.6;\n";
	ofile<<"camera {";
	frameCamera.GetOrigin(co);
	ofile<<"  location <"<<(-1)*(co[0])<<","<<co[1]<<","<<co[2]<<">";
	frameCamera.GetForwardVector(ca);
	ofile<<"  look_at <"<<(-1)*(co[0]+ca[0])<<","<<co[1]+ca[1]<<","<<co[2]+ca[2]<<">";
	ofile<<"  sky y";
	ofile<<"}\n";
	ofile<<"background{rgb 1}\n";
	ofile<<"\n";
	ofile<<"#declare REFLECTION = 0.02;\n#declare SPECULAR = 0.2;\n#declare AMBIENT = 0.78;\n#declare DIFFUSE = 0.8;\n";
	ofile<<"light_source {<0,0,1000>, color rgb 1 parallel point_at <0,0,0>}\n";

	for(int i=0; i<spheres.size(); i++)
	{
		M3DMatrix44f m0, m1, m2;
		AllSpheres.GetMatrix(m2);
		spheres[i].GetMatrix(m1);
		m3dMatrixMultiply44(m0, m2, m1);
		ofile<<"sphere {<"<<(-1)*(m0[12])<<","<<m0[13]<<","<<m0[14]<<">,"<<Spheres[i].radius<<" pigment { color rgb <";
		ofile<<Spheres[i].red<<","<<Spheres[i].green<<","<<Spheres[i].blue<<"> } finish{reflection REFLECTION specular SPECULAR ambient AMBIENT diffuse DIFFUSE} } \n";
	}
	for(int i=0; i<Lines.size(); i++)
	{
		M3DMatrix44f m0, m1, m2;
		AllSpheres.GetMatrix(m2);
		for(int j=0; j<12; j++) m1[j]=0.0;
		m1[12]=Lines[i].x1;
		m1[13]=Lines[i].y1;
		m1[14]=Lines[i].z1;
		m1[15]=0.0;
		m3dMatrixMultiply44(m0, m2, m1);
		ofile<<"cylinder {<"<<(-1)*(m0[12])<<","<<m0[13]<<","<<m0[14]<<">,";

		AllSpheres.GetMatrix(m2);
		for(int j=0; j<12; j++) m1[j]=0.0;
		m1[12]=Lines[i].x2;
		m1[13]=Lines[i].y2;
		m1[14]=Lines[i].z2;
		m1[15]=0.0;
		m3dMatrixMultiply44(m0, m2, m1);
		ofile<<"<"<<(-1)*(m0[12])<<","<<m0[13]<<","<<m0[14]<<">,"<<Spheres[0].radius*0.2<<" pigment { color rgb <";
		ofile<<Lines[i].red<<","<<Lines[i].green<<","<<Lines[i].blue<<"> } } \n";
	}
}
///////////////////////////////////////////////////////////////////////
// Draw spheres

//get depth of spheres[i]
/** \brief Get the camera-space depth of sphere i. (CLAUDE)
 *  \param i  Index of the sphere.
 *  \return   The depth (camera-space Z) of the sphere's center. */
double GetDepth(size_t i)
{
	M3DMatrix44f m0, m1, m2;
	frameCamera.GetCameraOrientation(m0);
	AllSpheres.GetMatrix(m1);
	m3dMatrixMultiply44(m2, m0, m1);
	spheres[i].GetMatrix(m1);
	m3dMatrixMultiply44(m0, m2, m1);
	return m0[14];
}
/** \brief Draw all spheres (depth-sorted) and lines of the current scene. (CLAUDE) */
void DrawInhabitants()
{
	// Draw the spheres
	{
		//sort the spheres according to the depth of their centers
		auto DepthCompareFunc = [&] (const size_t i, const size_t j) -> bool
		{
			return GetDepth(i)<GetDepth(j);
		};
		std::sort(RenderOrder.begin(), RenderOrder.end(), DepthCompareFunc);
	}
	glPushMatrix();
	AllSpheres.ApplyActorTransform();
	glEnable(GL_LIGHTING);
	for(auto iter=RenderOrder.begin(); iter!=RenderOrder.end(); iter++)
	{
		size_t i=*iter;
		glPushMatrix();
		glColor4f(Spheres.at(i).red, Spheres.at(i).green, Spheres.at(i).blue, Spheres.at(i).transparency);
		spheres[i].ApplyActorTransform();
		glutSolidSphere(Spheres.at(i).radius, Precision, Precision);
		glPopMatrix();
	}
	glDisable(GL_LIGHTING);
	for(std::vector<line>::const_iterator iter=Lines.begin(); iter!=Lines.end(); iter++)
	{
		glColor4f(iter->red, iter->green, iter->blue, iter->transparency);
		glLineWidth(iter->width);
		glBegin(GL_LINES);
		glVertex3d(iter->x1, iter->y1, iter->z1);
		glVertex3d(iter->x2, iter->y2, iter->z2);
		glEnd();
	}
	glPopMatrix();

}


// Called to draw scene
/** \brief GLUT display callback: render the 3D scene and 2D overlay, then swap buffers. */
void RenderScene(void)
{
	// Clear the window with current clearing color
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	frameCamera.ApplyCameraTransform();
	// Position light before any other transformations
	glLightfv(GL_LIGHT0, GL_POSITION, fLightPos);
	DrawInhabitants();
	glPopMatrix();


	//2d stuff
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, WindowWidth, 0, WindowHeight, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);

	(*pExtra2DStuffFunc)(WindowWidth, WindowHeight);

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	// Do the buffer Swap
	glutSwapBuffers();
	glutPostRedisplay();
}


static double OrthoWindowSize=1.0;
/** \brief Switch the projection to orthographic mode. (CLAUDE) */
void ChangeToOrtho()
{
	OrthoProjection=true;
	// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(fAspect*(-1)*OrthoWindowSize, fAspect*OrthoWindowSize, ((-1)*OrthoWindowSize), OrthoWindowSize, 0.1, 500);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();    
}
/** \brief Switch the projection to perspective mode. (CLAUDE) */
void ChangeToPerspective()
{
	OrthoProjection=false;
	// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(35.0f, fAspect, 0.1f, 500.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();    
}

// Respond to key strokes
/** \brief GLUT special-key callback: arrow keys move/zoom the camera. (CLAUDE)
 *  \param key  GLUT special key code.
 *  \param x    Mouse X position when the key was pressed.
 *  \param y    Mouse Y position when the key was pressed. */
void SpecialKeys(int key, int x, int y)
{
	if(key == GLUT_KEY_UP)
		if(OrthoProjection)
		{
			OrthoWindowSize*=0.9;
			ChangeToOrtho();
		}
		else
			frameCamera.MoveForward(0.1f);

	if(key == GLUT_KEY_DOWN)
		if(OrthoProjection)
		{
			OrthoWindowSize*=1.1;
			ChangeToOrtho();
		}
		else
			frameCamera.MoveForward(-0.1f);

	if(key == GLUT_KEY_LEFT)
		frameCamera.RotateLocalY(0.1);

	if(key == GLUT_KEY_RIGHT)
		frameCamera.RotateLocalY(-0.1);

	// Refresh the Window
	glutPostRedisplay();
}
/** \brief GLUT keyboard callback: handle precision, projection, movie, quit, and save keys. (CLAUDE)
 *  \param key  ASCII code of the pressed key.
 *  \param x    Mouse X position when the key was pressed.
 *  \param y    Mouse Y position when the key was pressed. */
void KeyStroke(unsigned char key, int x, int y)
{
	if(key == '+')
	{
		::Precision *= 2;
		if( ::Precision==0)
			::Precision=1;
		glutPostRedisplay();
	}
	if(key == '-')
	{
		::Precision /= 2;
		glutPostRedisplay();
	}
	if(key == 'q')
	{
		glutDestroyWindow( ::window_handle);
		glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);
		glutLeaveMainLoop();
	}
	if(key == 'o')
	{
		ChangeToOrtho();
		glutPostRedisplay();
	}
	if(key == 'p')
	{
		ChangeToPerspective();
		glutPostRedisplay();
	}
	if(key == 'a')
	{
		FrameRate/=1.5;
	}
	if(key == 's')
	{
		FrameRate*=1.5;
	}
	if(key == 'd')
	{
		Playing = !Playing;
		LastFrameTime = 0;
	}
	if(key == 'r')
	{
		Reverse = !Reverse;
	}
	if(key=='v')
		SavePovRay();
}

/** \brief GLUT idle callback: advance the movie frame when playing. (CLAUDE) */
void IdleFunction(void)
{
	if(Playing)
	{
		long Time = GetPreciseClock();
		if(Time-LastFrameTime>FrameRate || LastFrameTime==0)
		{
			int NextFrame = CurrentFrame;
			size_t FrameChange;
			if (LastFrameTime == 0)
				FrameChange = 1;
			else
				FrameChange = (Time - LastFrameTime) / FrameRate;
			if (Reverse)
				NextFrame -= FrameChange;
			else
				NextFrame += FrameChange;
			if(NextFrame>=0 && NextFrame<MaxFrame)
			{
				PrepareFrame(NextFrame);
				glutPostRedisplay();
			}
			else
			{
				Playing=false;
			}
		}
	}
}
///////////////////////////////////////////////////////////
/** \brief GLUT reshape callback: update viewport and projection for a new window size. (CLAUDE)
 *  \param w  New window width in pixels.
 *  \param h  New window height in pixels. */
void ChangeSize(int w, int h)
{

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if(h == 0)
		h = 1;

	glViewport(0, 0, w, h);

	fAspect = (GLfloat)w / (GLfloat)h;
	WindowWidth=w;
	WindowHeight=h;

	// Set the clipping volume
	if(OrthoProjection)
		ChangeToOrtho();
	else
		ChangeToPerspective();
}

///////////////////////////////////////////////////////////
//this part of code deals with mouse dragging
static int CurrentMouseX, CurrentMouseY;
/** \brief GLUT passive-motion callback: track the current mouse position. (CLAUDE)
 *  \param x  Mouse X position.
 *  \param y  Mouse Y position. */
void PassiveMotion(int x, int y)
{
	::CurrentMouseX=x;
	::CurrentMouseY=y;
}
float Sensitivity = 0.003;
/** \brief GLUT motion callback: rotate the scene by the mouse drag delta. (CLAUDE)
 *  \param x  Current mouse X position.
 *  \param y  Current mouse Y position. */
void Motion(int x, int y)
{
	int dx=x-CurrentMouseX;
	int dy=y-CurrentMouseY;
	::CurrentMouseX=x;
	::CurrentMouseY=y;

	AllSpheres.RotateWorld( ::Sensitivity*dy, 1.0, 0.0, 0.0);
	AllSpheres.RotateWorld( ::Sensitivity*dx, 0.0, 1.0, 0.0);
}

/** \brief Set up GLUT and run the main loop to display a multi-frame sphere/line scene. (CLAUDE)
 *  \param GetFrameFunc  Callback that fills spheres and lines for a given frame number.
 *  \param NumFrames     Total number of frames.
 *  \param Extra2DStuff  Callback receiving (window width, height) to draw extra 2D overlay content. */
void DisplaySpheres_MultipleFrames(std::function<void(std::vector<sphere> & VSpheres, std::vector<line> & VLines, long FrameNumber)> GetFrameFunc, long NumFrames, std::function<void(size_t, size_t)> Extra2DStuff)
{
	::MaxFrame=NumFrames;
	pGetFrameFunc = & GetFrameFunc;
	pExtra2DStuffFunc = & Extra2DStuff;


	int fake_argc=1;
	char * fake_argv="DisplaySpheres";
	glutInit(&fake_argc, &fake_argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(WindowWidth, WindowHeight);
	window_handle=glutCreateWindow("");
	glutReshapeFunc(ChangeSize);
	glutDisplayFunc(RenderScene);
	glutSpecialFunc(SpecialKeys);
	glutIdleFunc(IdleFunction);
	glutMotionFunc(Motion);
	glutPassiveMotionFunc(PassiveMotion);
	glutKeyboardFunc(KeyStroke);

	SetupRC();

	glutMainLoop();
}
void DisplaySpheres_Movie(std::function<void(std::vector<sphere> & VSpheres, std::vector<line> & VLines, long FrameNumber)> GetFrameFunc, long NumFrames)
{
	auto ExtraDisplayFunc = [&] (size_t, size_t) ->void
	{
		std::stringstream ss;
		ss<<"Frame "<<CurrentFrame;
		DisplayText(ss.str());
	};
	DisplaySpheres_MultipleFrames(GetFrameFunc, NumFrames, ExtraDisplayFunc);
}
void DisplaySpheres(const std::vector<sphere> & VSpheres, const std::vector<line> & VLines, std::function<void(size_t, size_t)> Extra2DStuff)
{
	auto GetFrameFunc = [&](std::vector<sphere> & V1, std::vector<line> & V2, long FrameNumber) ->void
	{
		V1=VSpheres;
		V2=VLines;
	};
	DisplaySpheres_MultipleFrames(GetFrameFunc, 1, Extra2DStuff);
}

//transferred from PeriodicCellList.h
/** \brief Embed a 1D, 2D, or 3D configuration into 3D space for display. (CLAUDE)
 *  \param a  Source configuration (dimension 1, 2, or 3).
 *  \return   An equivalent 3D configuration suitable for rendering. */
Configuration To3D(const Configuration & a)
{
	if (a.GetDimension() == 3)
	{
		return a;
	}
	else if (a.GetDimension() == 2)
	{
		GeometryVector bas[3];
		bas[0] = a.GetBasisVector(0);
		bas[0].SetDimension(3);
		bas[1] = a.GetBasisVector(1);
		bas[1].SetDimension(3);
		bas[2] = GeometryVector(0.0, 0.0, 0.1*std::sqrt(a.PeriodicVolume()));

		Configuration b(3, bas, 1.0);
		for (size_t i = 0; i<a.NumParticle(); i++)
		{
			//auto pa=a.GetParticle(i);
			GeometryVector r = a.GetRelativeCoordinates(i);
			r.SetDimension(3);
			r.x[2] = 0.5;
			b.Insert(a.GetCharacteristics(i), r);
	}

		return b;
}
	else if (a.GetDimension() == 1)
	{
		GeometryVector bas[3];
		bas[0] = a.GetBasisVector(0);
		bas[0].SetDimension(3);
		bas[1] = GeometryVector(0.0, 0.1*a.PeriodicVolume(), 0.0);
		bas[2] = GeometryVector(0.0, 0.0, 0.1*a.PeriodicVolume());

		Configuration b(3, bas, 1.0);
		for (size_t i = 0; i<a.NumParticle(); i++)
		{
			//auto pa=a.GetParticle(i);
			GeometryVector r = a.GetRelativeCoordinates(i);
			r.SetDimension(3);
			r.x[2] = 0.5;
			r.x[1] = 0.5;
			b.Insert(a.GetCharacteristics(i), r);
		}

		return b;
	}
	else
	{
		std::cerr << "Error in PeriodicCellList.cpp : To3D() : unsupported dimension!\n";
		return GetUnitCubicBox(3);
	}
}
/** \brief Convert a configuration into renderable spheres (particles), bonds, and box edges. (CLAUDE)
 *  \param list           Source 3D configuration.
 *  \param spheres        Destination vector receiving the particle spheres.
 *  \param lines          Destination vector receiving box-edge and bond lines.
 *  \param Radius         Sphere radius; 0 selects a radius from typicalLength.
 *  \param BondLength     Maximum neighbor distance to draw bonds; <= 0 draws no bonds.
 *  \param typicalLength  Typical inter-particle spacing used to auto-size spheres.
 *  \param OriginalDim    Original dimensionality of the configuration before embedding in 3D. */
void ConfigurationtoSpheres(const Configuration & list, std::vector<sphere> & spheres, std::vector<line> & lines, double Radius, double BondLength, double typicalLength, DimensionType OriginalDim)
{
	std::vector<GeometryVector> a(8, GeometryVector(0.0, 0.0, 0.0));
	a[1] = a[0] + list.GetBasisVector(0);
	a[2] = a[0] + list.GetBasisVector(1);
	a[3] = a[0] + list.GetBasisVector(2);
	a[4] = a[1] + list.GetBasisVector(1);
	a[5] = a[1] + list.GetBasisVector(2);
	a[6] = a[2] + list.GetBasisVector(2);
	a[7] = a[4] + list.GetBasisVector(2);
	GeometryVector DisplayShift = a[7] * 0.5;
	double one_over_scale = 1.0 / std::sqrt(DisplayShift.Modulus2());
	for (int i = 0; i<8; i++)
	{
		a[i].MinusFrom(DisplayShift);
		a[i].MultiplyFrom(one_over_scale);
	}
	int point1[12] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 5, 6 };
	int point2[12] = { 1, 2, 3, 4, 5, 4, 6, 5, 6, 7, 7, 7 };
	for (int i = 0; i<12; i++)
	{
		line temp;
		temp.x1 = a[point1[i]].x[0];
		temp.y1 = a[point1[i]].x[1];
		temp.z1 = a[point1[i]].x[2];
		temp.x2 = a[point2[i]].x[0];
		temp.y2 = a[point2[i]].x[1];
		temp.z2 = a[point2[i]].x[2];
		temp.blue = 0.0;
		temp.red = 0.0;
		temp.green = 0.0;
		temp.transparency = 1.0;
		temp.width = 3;
		lines.push_back(temp);
	}

	for (size_t i = 0; i < list.NumParticle(); i++)
	{
		GeometryVector loc = (list.GetCartesianCoordinates(i) - DisplayShift)*one_over_scale;
		sphere temp;
		temp.x = loc.x[0];
		temp.y = loc.x[1];
		temp.z = loc.x[2];

		char chr1 = list.GetCharacteristics(i).name[0];
		if (chr1 >= 'A' && chr1 <= 'Z')
		{
			temp.red = ((chr1 - 'A') / 6)*0.1 + 0.3;
			temp.green = ((chr1 - 'A') % 6)*0.1 + 0.3;
			char chr2 = list.GetCharacteristics(i).name[1];
			if (chr2 >= 'a' && chr2 <= 'z')
				temp.blue = (chr2 - 'a')*0.02 + 0.4;
			else
				temp.blue = 0.35;
			//debug temp
			//if(chr1=='A')
			//{
			//	temp.blue=1.0;
			//	temp.red=1.0;
			//	temp.green=1.0;
			//}
		}
		else
		{
			temp.red = 0.25;
			temp.green = 0.25;
			temp.blue = 0.25;
		}
		temp.transparency = 1.00;

		if (Radius == 0.0)
			temp.radius = 0.3*typicalLength*one_over_scale;
		else
			temp.radius = Radius * one_over_scale;

		spheres.push_back(temp);

		if (BondLength > 0)
		{
			//debug temp
			size_t a = lines.size();

			auto IterateFunc = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
			{
				if (shift.Modulus2() == 0.0)
					return;
				for (DimensionType i = OriginalDim; i < 3; i++)
					if (PeriodicShift[i] != 0)
						return;
				GeometryVector end = loc + shift * one_over_scale;
				line temp;
				temp.x1 = loc.x[0];
				temp.y1 = loc.x[1];
				temp.z1 = loc.x[2];
				temp.x2 = end.x[0];
				temp.y2 = end.x[1];
				temp.z2 = end.x[2];
				temp.blue = 0.25;
				temp.red = 0.25;
				temp.green = 0.25;
				temp.transparency = 1.0;
				temp.width = 2.0;
				lines.push_back(temp);
			};
			list.IterateThroughNeighbors(list.GetRelativeCoordinates(i), BondLength, IterateFunc);

			//debug temp
			//std::cout << lines.size() - a << '\n';
		}
	}

}
/** \brief Estimate a typical inter-particle spacing for a configuration. (CLAUDE)
 *  \param list  Configuration to measure.
 *  \return      A typical length scale (nearest-neighbor distance, bounded by density). */
double findTypicalLength(const Configuration & list)
{
	double typicalLength = std::pow(list.PeriodicVolume() / list.NumParticle(), 1.0 / list.GetDimension());
	double MinLength = 3 * typicalLength;//distance between nearest particles
	for (size_t i = 0; i < list.NumParticle(); i++)
	{
		//Configuration::particle * pa=list.GetParticle(i);
		list.IterateThroughNeighbors(list.GetRelativeCoordinates(i), MinLength, [&MinLength, &i](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
		{
			if (SourceAtom != i)
			{
				double l = std::sqrt(shift.Modulus2());
				if (l < MinLength)
					MinLength = l;
			}
		});
	}
	if (MinLength > 0.3*typicalLength && MinLength < typicalLength)
		//we could make sure no spheres overlap (desirable!) by slightly decreasing sphere radius, 
		//so do it!
		typicalLength = MinLength;
	return typicalLength;
}

void DisplayConfiguration(const Configuration & a, double BondLength, double Radius)
{
	std::vector<sphere> spheres;
	std::vector<line> lines;

	ConfigurationtoSpheres(To3D(a), spheres, lines, Radius, BondLength, findTypicalLength(a), a.GetDimension());


	DisplaySpheres(spheres, lines);
}

void Display3DConfigurationMovie(ConfigurationPack & pk, double BondLength, double Radius)
{
	if (pk.NumConfig() == 0)
	{
		std::cout << "warning in Display3DConfigurationMovie : ConfigurationPack is empty. Returning!\n";
		return;
	}
	Configuration a = pk.GetConfig(0);

	double typicalLength = findTypicalLength(a);
	auto GetFrameFunc = [&a, &pk, &Radius, &BondLength, &typicalLength](std::vector<sphere> & spheres, std::vector<line> & lines, long FrameNumber) ->void
	{
		spheres.clear();
		lines.clear();
		Configuration l = To3D(pk.GetConfig(FrameNumber));
		ConfigurationtoSpheres(l, spheres, lines, Radius, BondLength, typicalLength, a.GetDimension());
	};

	DisplaySpheres_Movie(GetFrameFunc, pk.NumConfig());
}

void Display3DStructureFactor(const Configuration& a, double KCutOff, double KSphereRadius)
{
	assert(a.GetDimension() == 3);
	double Smin = 1000000, Smax = -1;
	auto SetColorFunc = [&](double Sk, sphere & s) ->void
	{
		s.red = (Sk - Smin) / (Smax - Smin);
		s.green = 0.0;
		s.blue = 1.0 - s.red;
	};
	auto DrawColorBarFunc = [&](size_t w, size_t h) ->void
	{
		glBegin(GL_POLYGON);
		glColor3f(0.0, 0.0, 1.0);
		glVertex2f(0.87*w, 0.1*h);
		glVertex2f(0.9*w, 0.1*h);
		glColor3f(1.0, 0.0, 0.0);
		glVertex2f(0.9*w, 0.9*h);
		glVertex2f(0.87*w, 0.9*h);
		glEnd();
		std::stringstream smin;
		smin << Smin;
		DisplayText(smin.str(), 0.91*w, 0.1*h);
		std::stringstream smax;
		smax << Smax;
		DisplayText(smax.str(), 0.91*w, 0.9*h);
};
	GeometryVector recip[3];
	for (int i = 0; i<3; i++)
		recip[i] = a.GetReciprocalBasisVector(i);
	PeriodicCellList<Empty> lReciprocal(3, recip, MaxDistance, false);
	lReciprocal.Insert(Empty(), GeometryVector(3));
	std::vector<sphere> VSpheres;
	auto IterateFunc = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle)
	{
		if (shift.Modulus2() == 0.0)
			return;
		sphere s;
		s.x = shift.x[0] / KCutOff;
		s.y = shift.x[1] / KCutOff;
		s.z = shift.x[2] / KCutOff;
		double sk = std::log10(StructureFactor(a, shift));
		if (sk<Smin)
			Smin = sk;
		if (sk>Smax)
			Smax = sk;
		s.transparency = sk;//temperarily use this data member to store S(k)
		VSpheres.push_back(s);
	};
	lReciprocal.IterateThroughNeighbors(GeometryVector(3), KCutOff, IterateFunc);
	if (VSpheres.empty())
		return;
	Smax = std::ceil(Smax);
	Smin = std::floor(Smin);
	if (KSphereRadius == 0.0)
	{
		//determine the correct radius of spheres
		auto CompareFunc = [](const sphere & left, const sphere & right) ->bool
		{
			return left.x*left.x + left.y*left.y + left.z*left.z<right.x*right.x + right.y*right.y + right.z*right.z;
		};
		std::nth_element(VSpheres.begin(), VSpheres.begin(), VSpheres.end(), CompareFunc);
		sphere s = VSpheres[0];
		KSphereRadius = 0.2 * std::sqrt(s.x*s.x + s.y*s.y + s.z*s.z);
	}
	for (auto iter = VSpheres.begin(); iter != VSpheres.end(); iter++)
	{
		SetColorFunc(iter->transparency, *iter);

		//set transparency accoring to intensity
		//this enables us to view non-zero values better
		iter->transparency = (iter->transparency - Smin) / (Smax - Smin)*0.8 + 0.2;
		//iter->transparency=0.5;

		iter->radius = KSphereRadius;
	}
	{
		//use a black sphere to indicate origin
		sphere s = VSpheres[0];
		s.x = s.y = s.z = s.red = s.green = s.blue = 0.0;
		s.transparency = 1.0;
		VSpheres.push_back(s);
	}
	DisplaySpheres(VSpheres, std::vector<line>(), DrawColorBarFunc);
}





#else

//These functions are disabled if OpenGL is disabled
void DisplayText(const std::string & text, double x, double y)
{
	return;
}
void DisplaySpheres(const std::vector<sphere> & VSpheres, const std::vector<line> & VLines, std::function<void(size_t, size_t)> Extra2DStuff)
{
	return;
}
void DisplaySpheres_Movie(std::function<void(std::vector<sphere> & VSpheres, std::vector<line> & VLines, long FrameNumber)> GetFrameFunc, long NbrFrames)
{
	return;
}

void Display3DStructureFactor(const Configuration& a, double KCutOff, double KSphereRadius)
{
	std::cerr << "Error in Display3DStructureFactor : OpenGL not enabled!\n";
}

void Display3DConfigurationMovie(ConfigurationPack & pk, double BondLength, double Radius) 
{
	std::cerr << "Error in Display3DConfigurationMovie : OpenGL not enabled!\n";
}

void DisplayConfiguration(const Configuration & a, double BondLength, double Radius) 
{
	std::cerr << "Error in DisplayConfiguration : OpenGL not enabled!\n";
}


#endif /* USE_OPENGL*/
