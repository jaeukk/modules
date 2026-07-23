// Frame.h
// Implementation of the GLFrame Class
// Richard S. Wright Jr.
// Code by Richard S. Wright Jr.

#include "math3d.h"

#ifdef USE_OPENGL

#include <GL/freeglut.h>

#ifndef _ORTHO_FRAME_
#define _ORTHO_FRAME_

// The GLFrame (OrthonormalFrame) class. Possibly the most useful little piece of 3D graphics
// code for OpenGL immersive environments.
// Richard S. Wright Jr.
class GLFrame
    {
	protected:
        M3DVector3f vOrigin;	// Where am I?
        M3DVector3f vForward;	// Where am I going?
        M3DVector3f vUp;		// Which way is up?

    public:
		// Default position and orientation. At the origin, looking
		// down the positive Z axis (right handed coordinate system).
		/** \brief Construct a frame at the origin looking down -Z with +Y up. */
		GLFrame(void) {
			// At origin
            vOrigin[0] = 0.0f; vOrigin[1] = 0.0f; vOrigin[2] = 0.0f; 

			// Up is up (+Y)
            vUp[0] = 0.0f; vUp[1] = 1.0f; vUp[2] = 0.0f;

			// Forward is -Z (default OpenGL)
            vForward[0] = 0.0f; vForward[1] = 0.0f; vForward[2] = -1.0f;
            }


        /////////////////////////////////////////////////////////////
        // Set Location
        /** \brief Set the frame origin from a 3D point. (CLAUDE)
         *  \param vPoint  New origin location. */
        inline void SetOrigin(const M3DVector3f vPoint) {
			m3dCopyVector3(vOrigin, vPoint); }

        /** \brief Set the frame origin from individual coordinates. (CLAUDE)
         *  \param x  X coordinate of the origin.
         *  \param y  Y coordinate of the origin.
         *  \param z  Z coordinate of the origin. */
        inline void SetOrigin(float x, float y, float z) {
			vOrigin[0] = x; vOrigin[1] = y; vOrigin[2] = z; }

		/** \brief Get the frame origin into a 3D vector. (CLAUDE)
		 *  \param vPoint  Destination 3D vector receiving the origin. */
		inline void GetOrigin(M3DVector3f vPoint) {
			m3dCopyVector3(vPoint, vOrigin); }

		/** \brief Return the X coordinate of the frame origin. (CLAUDE) */
		inline float GetOriginX(void) { return vOrigin[0]; }
		/** \brief Return the Y coordinate of the frame origin. (CLAUDE) */
		inline float GetOriginY(void) { return vOrigin[1]; }
		/** \brief Return the Z coordinate of the frame origin. (CLAUDE) */
		inline float GetOriginZ(void) { return vOrigin[2]; }

        /////////////////////////////////////////////////////////////
        // Set Forward Direction
        /** \brief Set the forward (view) direction from a 3D vector. (CLAUDE)
         *  \param vDirection  New forward direction. */
        inline void SetForwardVector(const M3DVector3f vDirection) {
			m3dCopyVector3(vForward, vDirection); }

        /** \brief Set the forward (view) direction from individual components. (CLAUDE)
         *  \param x  X component of the forward direction.
         *  \param y  Y component of the forward direction.
         *  \param z  Z component of the forward direction. */
        inline void SetForwardVector(float x, float y, float z)
            { vForward[0] = x; vForward[1] = y; vForward[2] = z; }

        /** \brief Get the forward (view) direction into a 3D vector. (CLAUDE)
         *  \param vVector  Destination 3D vector receiving the forward direction. */
        inline void GetForwardVector(M3DVector3f vVector) { m3dCopyVector3(vVector, vForward); }

        /////////////////////////////////////////////////////////////
        // Set Up Direction
        /** \brief Set the up direction from a 3D vector. (CLAUDE)
         *  \param vDirection  New up direction. */
        inline void SetUpVector(const M3DVector3f vDirection) {
			m3dCopyVector3(vUp, vDirection); }

        /** \brief Set the up direction from individual components. (CLAUDE)
         *  \param x  X component of the up direction.
         *  \param y  Y component of the up direction.
         *  \param z  Z component of the up direction. */
        inline void SetUpVector(float x, float y, float z)
			{ vUp[0] = x; vUp[1] = y; vUp[2] = z; }

        /** \brief Get the up direction into a 3D vector. (CLAUDE)
         *  \param vVector  Destination 3D vector receiving the up direction. */
        inline void GetUpVector(M3DVector3f vVector) { m3dCopyVector3(vVector, vUp); }


		/////////////////////////////////////////////////////////////
		// Get Axes
		/** \brief Get the local Z axis (the forward vector) into a 3D vector. (CLAUDE)
		 *  \param vVector  Destination 3D vector receiving the Z axis. */
		inline void GetZAxis(M3DVector3f vVector) { GetForwardVector(vVector); }
		/** \brief Get the local Y axis (the up vector) into a 3D vector. (CLAUDE)
		 *  \param vVector  Destination 3D vector receiving the Y axis. */
		inline void GetYAxis(M3DVector3f vVector) { GetUpVector(vVector); }
		/** \brief Get the local X axis (up cross forward) into a 3D vector. (CLAUDE)
		 *  \param vVector  Destination 3D vector receiving the X axis. */
		inline void GetXAxis(M3DVector3f vVector) { m3dCrossProduct(vVector, vUp, vForward); }


		/////////////////////////////////////////////////////////////
        // Translate along orthonormal axis... world or local
        /** \brief Translate the origin in world coordinates. (CLAUDE)
         *  \param x  Offset along world X.
         *  \param y  Offset along world Y.
         *  \param z  Offset along world Z. */
        inline void TranslateWorld(float x, float y, float z)
			{ vOrigin[0] += x; vOrigin[1] += y; vOrigin[2] += z; }

        /** \brief Translate the origin in local (frame) coordinates. (CLAUDE)
         *  \param x  Offset along the local right (X) axis.
         *  \param y  Offset along the local up (Y) axis.
         *  \param z  Offset along the local forward (Z) axis. */
        inline void TranslateLocal(float x, float y, float z)
			{ MoveForward(z); MoveUp(y); MoveRight(x);	}


		/////////////////////////////////////////////////////////////
		// Move Forward (along Z axis)
		/** \brief Move the origin along the forward (Z) axis. (CLAUDE)
		 *  \param fDelta  Distance to move along the forward direction. */
		inline void MoveForward(float fDelta)
			{
		    // Move along direction of front direction
			vOrigin[0] += vForward[0] * fDelta;
			vOrigin[1] += vForward[1] * fDelta;
			vOrigin[2] += vForward[2] * fDelta;
			}

		// Move along Y axis
		/** \brief Move the origin along the up (Y) axis. (CLAUDE)
		 *  \param fDelta  Distance to move along the up direction. */
		inline void MoveUp(float fDelta)
			{
		    // Move along direction of up direction
			vOrigin[0] += vUp[0] * fDelta;
			vOrigin[1] += vUp[1] * fDelta;
			vOrigin[2] += vUp[2] * fDelta;
			}

		// Move along X axis
		/** \brief Move the origin along the right (X) axis. (CLAUDE)
		 *  \param fDelta  Distance to move along the right direction. */
		inline void MoveRight(float fDelta)
			{
			// Move along direction of right vector
			M3DVector3f vCross;
			m3dCrossProduct(vCross, vUp, vForward);

			vOrigin[0] += vCross[0] * fDelta;
			vOrigin[1] += vCross[1] * fDelta;
			vOrigin[2] += vCross[2] * fDelta;
			}

		///////////////////////////////////////////////////////////////////////
		// Just assemble the 






		/** \brief Assemble the 4x4 transformation matrix for this frame. (CLAUDE)
		 *  \param matrix        Destination 4x4 matrix.
		 *  \param bRotationOnly  If true, omit the translation (rotation only). */
		void GetMatrix(M3DMatrix44f	matrix, bool bRotationOnly = false)
			{
			// Calculate the right side (x) vector, drop it right into the matrix
			M3DVector3f vXAxis;
			m3dCrossProduct(vXAxis, vUp, vForward);

			// Set matrix column does not fill in the fourth value...
            m3dSetMatrixColumn44(matrix, vXAxis, 0);
            matrix[3] = 0.0f;
           
            // Y Column
			m3dSetMatrixColumn44(matrix, vUp, 1);
            matrix[7] = 0.0f;       
                                    
            // Z Column
			m3dSetMatrixColumn44(matrix, vForward, 2);
            matrix[11] = 0.0f;

            // Translation (already done)
			if(bRotationOnly == true)
				{
				matrix[12] = 0.0f;
				matrix[13] = 0.0f;
				matrix[14] = 0.0f;
				}
			else
				m3dSetMatrixColumn44(matrix, vOrigin, 3);

            matrix[15] = 1.0f;
			}


        /////////////////////////////////////////////////////////////
        // Get a 4x4 transformation matrix that describes the ccamera
        // orientation.
        /** \brief Build the 4x4 camera orientation matrix (rotation only, no translation).
         *  \param m  Destination 4x4 camera orientation matrix. */
        inline void GetCameraOrientation(M3DMatrix44f m)
            {
			M3DVector3f x, z;

			// Make rotation matrix
			// Z vector is reversed
			z[0] = -vForward[0];
			z[1] = -vForward[1];
			z[2] = -vForward[2];

			// X vector = Y cross Z 
			m3dCrossProduct(x, vUp, z);

			// Matrix has no translation information and is
			// transposed.... (rows instead of columns)
			#define M(row,col)  m[col*4+row]
			   M(0, 0) = x[0];
			   M(0, 1) = x[1];
			   M(0, 2) = x[2];
			   M(0, 3) = 0.0;
			   M(1, 0) = vUp[0];
			   M(1, 1) = vUp[1];
			   M(1, 2) = vUp[2];
			   M(1, 3) = 0.0;
			   M(2, 0) = z[0];
			   M(2, 1) = z[1];
			   M(2, 2) = z[2];
			   M(2, 3) = 0.0;
			   M(3, 0) = 0.0;
			   M(3, 1) = 0.0;
			   M(3, 2) = 0.0;
			   M(3, 3) = 1.0;
			#undef M 
            }
            
            
		/////////////////////////////////////////////////////////////
		// Perform viewing or modeling transformations
		// Position as the camera (for viewing). Apply this transformation
		// first as your viewing transformation
		// The default implementation of gluLookAt can be considerably sped up
		// since it uses doubles for everything... then again profile before you
		// tune... ;-) You might get a boost form page fault reduction too... if
		// no other glu routines are used...
		// This will get called once per frame.... go ahead and inline
        /** \brief Apply this frame as the camera (viewing) transform to the GL matrix stack. (CLAUDE)
         *  \param bRotOnly  If true, apply rotation only and skip the translation. */
        inline void ApplyCameraTransform(bool bRotOnly = false)
			{
			M3DMatrix44f m;

            GetCameraOrientation(m);
            
			// Camera Transform   
			glMultMatrixf(m);
		
			// If Rotation only, then do not do the translation
			if(!bRotOnly)
				glTranslatef(-vOrigin[0], -vOrigin[1], -vOrigin[2]);

			/*gluLookAt(vOrigin[0], vOrigin[1], vOrigin[2],
						vOrigin[0] + vForward[0], 
						vOrigin[1] + vForward[1], 
						vOrigin[2] + vForward[2], 
						vUp[0], vUp[1], vUp[2]);
			*/
			}


		// Position as an object in the scene. This places and orients a
		// coordinate frame for other objects (besides the camera)
		// There is ample room for optimization here... 
		// This is going to be called alot... don't inline
		// Add flag to perform actor rotation only and not the translation
        /** \brief Apply this frame as an object (actor) transform to the GL matrix stack. (CLAUDE)
         *  \param bRotationOnly  If true, apply rotation only and skip the translation. */
        void ApplyActorTransform(bool bRotationOnly = false)
			{
			M3DMatrix44f rotMat;
			
			GetMatrix(rotMat, bRotationOnly);

			// Apply rotation to the current matrix
			glMultMatrixf(rotMat);
			}


        // Rotate around local X Axes - Note all rotations are in radians
        /** \brief Rotate the frame about its local X axis. (CLAUDE)
         *  \param fAngle  Rotation angle in radians. */
        void RotateLocalX(float fAngle)
			{
			M3DMatrix44f rotMat;
			M3DVector3f vCross;
			m3dCrossProduct(vCross, vUp, vForward);
			m3dRotationMatrix44(rotMat, fAngle,
							vCross[0], vCross[1], vCross[2]);

			M3DVector3f newVect;
			// Inline 3x3 matrix multiply for rotation only
			newVect[0] = rotMat[0] * vForward[0] + rotMat[4] * vForward[1] + rotMat[8] *  vForward[2];	
			newVect[1] = rotMat[1] * vForward[0] + rotMat[5] * vForward[1] + rotMat[9] *  vForward[2];	
			newVect[2] = rotMat[2] * vForward[0] + rotMat[6] * vForward[1] + rotMat[10] * vForward[2];	
			m3dCopyVector3(vForward, newVect);

			// Update pointing up vector
			newVect[0] = rotMat[0] * vUp[0] + rotMat[4] * vUp[1] + rotMat[8] *  vUp[2];	
			newVect[1] = rotMat[1] * vUp[0] + rotMat[5] * vUp[1] + rotMat[9] *  vUp[2];	
			newVect[2] = rotMat[2] * vUp[0] + rotMat[6] * vUp[1] + rotMat[10] * vUp[2];	
			m3dCopyVector3(vUp, newVect);
			}

		// Rotate around local Y
		/** \brief Rotate the frame about its local Y axis. (CLAUDE)
		 *  \param fAngle  Rotation angle in radians. */
        void RotateLocalY(float fAngle)
			{
	        M3DMatrix44f rotMat;

			// Just Rotate around the up vector
			// Create a rotation matrix around my Up (Y) vector
			m3dRotationMatrix44(rotMat, fAngle,
                         vUp[0], vUp[1], vUp[2]);

			M3DVector3f newVect;

	        // Rotate forward pointing vector (inlined 3x3 transform)
			newVect[0] = rotMat[0] * vForward[0] + rotMat[4] * vForward[1] + rotMat[8] *  vForward[2];	
			newVect[1] = rotMat[1] * vForward[0] + rotMat[5] * vForward[1] + rotMat[9] *  vForward[2];	
			newVect[2] = rotMat[2] * vForward[0] + rotMat[6] * vForward[1] + rotMat[10] * vForward[2];	
			m3dCopyVector3(vForward, newVect);
			}


		// Rotate around local Z
		/** \brief Rotate the frame about its local Z axis. (CLAUDE)
		 *  \param fAngle  Rotation angle in radians. */
        void RotateLocalZ(float fAngle)
			{
			M3DMatrix44f rotMat;

			// Only the up vector needs to be rotated
			m3dRotationMatrix44(rotMat, fAngle,
							vForward[0], vForward[1], vForward[2]);

			M3DVector3f newVect;
			newVect[0] = rotMat[0] * vUp[0] + rotMat[4] * vUp[1] + rotMat[8] *  vUp[2];	
			newVect[1] = rotMat[1] * vUp[0] + rotMat[5] * vUp[1] + rotMat[9] *  vUp[2];	
			newVect[2] = rotMat[2] * vUp[0] + rotMat[6] * vUp[1] + rotMat[10] * vUp[2];	
			m3dCopyVector3(vUp, newVect);
			}


		// Reset axes to make sure they are orthonormal. This should be called on occasion
		// if the matrix is long-lived and frequently transformed.
		/** \brief Re-orthonormalize the frame's up and forward axes. (CLAUDE) */
		void Normalize(void)
			{
			M3DVector3f vCross;

			// Calculate cross product of up and forward vectors
			m3dCrossProduct(vCross, vUp, vForward);

			// Use result to recalculate forward vector
			m3dCrossProduct(vForward, vCross, vUp);	

			// Also check for unit length...
			m3dNormalizeVector(vUp);
			m3dNormalizeVector(vForward);
			}


		// Rotate in world coordinates...
		/** \brief Rotate the frame about an arbitrary world-space axis. (CLAUDE)
		 *  \param fAngle  Rotation angle in radians.
		 *  \param x       X component of the world-space rotation axis.
		 *  \param y       Y component of the world-space rotation axis.
		 *  \param z       Z component of the world-space rotation axis. */
		void RotateWorld(float fAngle, float x, float y, float z)
			{
            M3DMatrix44f rotMat;

			// Create the Rotation matrix
			m3dRotationMatrix44(rotMat, fAngle, x, y, z);

			M3DVector3f newVect;
			
			// Transform the up axis (inlined 3x3 rotation)
			newVect[0] = rotMat[0] * vUp[0] + rotMat[4] * vUp[1] + rotMat[8] *  vUp[2];	
			newVect[1] = rotMat[1] * vUp[0] + rotMat[5] * vUp[1] + rotMat[9] *  vUp[2];	
			newVect[2] = rotMat[2] * vUp[0] + rotMat[6] * vUp[1] + rotMat[10] * vUp[2];	
			m3dCopyVector3(vUp, newVect);

			// Transform the forward axis
			newVect[0] = rotMat[0] * vForward[0] + rotMat[4] * vForward[1] + rotMat[8] *  vForward[2];	
			newVect[1] = rotMat[1] * vForward[0] + rotMat[5] * vForward[1] + rotMat[9] *  vForward[2];	
			newVect[2] = rotMat[2] * vForward[0] + rotMat[6] * vForward[1] + rotMat[10] * vForward[2];	
			m3dCopyVector3(vForward, newVect);
            }


        // Rotate around a local axis
        /** \brief Rotate the frame about an arbitrary local-space axis. (CLAUDE)
         *  \param fAngle  Rotation angle in radians.
         *  \param x       X component of the local-space rotation axis.
         *  \param y       Y component of the local-space rotation axis.
         *  \param z       Z component of the local-space rotation axis. */
        void RotateLocal(float fAngle, float x, float y, float z)
            {
            M3DVector3f vWorldVect;
			M3DVector3f vLocalVect;
			m3dLoadVector3(vLocalVect, x, y, z);

            LocalToWorld(vLocalVect, vWorldVect);
            RotateWorld(fAngle, vWorldVect[0], vWorldVect[1], vWorldVect[2]);
            }
    

		// Convert Coordinate Systems
        // This is pretty much, do the transformation represented by the rotation
        // and position on the point
		// Is it better to stick to the convention that the destination always comes
		// first, or use the conventions that "sounds" like the function...
        /** \brief Convert a point from local (frame) coordinates to world coordinates. (CLAUDE)
         *  \param vLocal  Source point in local coordinates.
         *  \param vWorld  Destination 3D vector receiving the point in world coordinates. */
        void LocalToWorld(const M3DVector3f vLocal, M3DVector3f vWorld)
            {
             // Create the rotation matrix based on the vectors
			M3DMatrix44f rotMat;

			GetMatrix(rotMat, true);

			// Do the rotation (inline it, and remove 4th column...)
			vWorld[0] = rotMat[0] * vLocal[0] + rotMat[4] * vLocal[1] + rotMat[8] *  vLocal[2];	
			vWorld[1] = rotMat[1] * vLocal[0] + rotMat[5] * vLocal[1] + rotMat[9] *  vLocal[2];	
			vWorld[2] = rotMat[2] * vLocal[0] + rotMat[6] * vLocal[1] + rotMat[10] * vLocal[2];	

            // Translate the point
            vWorld[0] += vOrigin[0];
            vWorld[1] += vOrigin[1];
            vWorld[2] += vOrigin[2];
            }

		// Change world coordinates into "local" coordinates
        /** \brief Convert a point from world coordinates to local (frame) coordinates. (CLAUDE)
         *  \param vWorld  Source point in world coordinates.
         *  \param vLocal  Destination 3D vector receiving the point in local coordinates. */
        void WorldToLocal(const M3DVector3f vWorld, M3DVector3f vLocal)
            {
			////////////////////////////////////////////////
            // Translate the origin
			M3DVector3f vNewWorld;
            vNewWorld[0] = vWorld[0] - vOrigin[0];
            vNewWorld[1] = vWorld[1] - vOrigin[1];
            vNewWorld[2] = vWorld[2] - vOrigin[2];

            // Create the rotation matrix based on the vectors
			M3DMatrix44f rotMat;
            M3DMatrix44f invMat;
			GetMatrix(rotMat, true);

			// Do the rotation based on inverted matrix
            m3dInvertMatrix44(invMat, rotMat);

			vLocal[0] = invMat[0] * vNewWorld[0] + invMat[4] * vNewWorld[1] + invMat[8] *  vNewWorld[2];	
			vLocal[1] = invMat[1] * vNewWorld[0] + invMat[5] * vNewWorld[1] + invMat[9] *  vNewWorld[2];	
			vLocal[2] = invMat[2] * vNewWorld[0] + invMat[6] * vNewWorld[1] + invMat[10] * vNewWorld[2];	
            }
        
        /////////////////////////////////////////////////////////////////////////////
        // Transform a point by frame matrix
        /** \brief Transform a point by the frame's full matrix (rotation and translation). (CLAUDE)
         *  \param vPointSrc  Source point.
         *  \param vPointDst  Destination 3D vector receiving the transformed point. */
        void TransformPoint(M3DVector3f vPointSrc, M3DVector3f vPointDst)
            {
            M3DMatrix44f m;
            GetMatrix(m, false);    // Rotate and translate
            vPointDst[0] = m[0] * vPointSrc[0] + m[4] * vPointSrc[1] + m[8] *  vPointSrc[2] + m[12];// * v[3];	 
            vPointDst[1] = m[1] * vPointSrc[0] + m[5] * vPointSrc[1] + m[9] *  vPointSrc[2] + m[13];// * v[3];	
            vPointDst[2] = m[2] * vPointSrc[0] + m[6] * vPointSrc[1] + m[10] * vPointSrc[2] + m[14];// * v[3];	
            }
        
        ////////////////////////////////////////////////////////////////////////////
        // Rotate a vector by frame matrix
        /** \brief Rotate a vector by the frame's matrix (rotation only, no translation). (CLAUDE)
         *  \param vVectorSrc  Source vector.
         *  \param vVectorDst  Destination 3D vector receiving the rotated vector. */
        void RotateVector(M3DVector3f vVectorSrc, M3DVector3f vVectorDst)
            {
            M3DMatrix44f m;
            GetMatrix(m, true);    // Rotate only
            
            vVectorDst[0] = m[0] * vVectorSrc[0] + m[4] * vVectorSrc[1] + m[8] *  vVectorSrc[2];	 
            vVectorDst[1] = m[1] * vVectorSrc[0] + m[5] * vVectorSrc[1] + m[9] *  vVectorSrc[2];	
            vVectorDst[2] = m[2] * vVectorSrc[0] + m[6] * vVectorSrc[1] + m[10] * vVectorSrc[2];	
            }
        };


#endif
#endif