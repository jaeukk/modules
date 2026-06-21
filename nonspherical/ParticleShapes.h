/**
*	Author	: Jaeuk Kim
*	Email	: phy000.kim@gmail.com	
*	Date	:	Feb. 2019 */

/** \file ParticleShapes.h
 * \brief Header file for the Shape and its child classes.
 *	An abstract struct "Shape" and its derivative structs are defined.
 *	These structs are mainly defined for computing the form factors.
 * 	However, they can be extented to other purposes, e.g., RSA procedure.
 *
 *	Shape --
 *			|-> PlanarFigure: Polygonal face in 3D
 *			|	|						|
 *			|-----> Polygon: 2D objects	|
 *			|							|
 *			|------------------------------> Polyhedron: 3D objects */


#ifndef PARTICLESHAPE_INCLUDED
#define PARTICLESHAPE_INCLUDED

#include <GeometryVector.h>
#include <RandomGenerator.h>
#include <etc.h>

#include <vector>
#include <complex>
#include <memory>
#include <bits/stdc++.h>

 /* To do list
 2. member function: check overlapping of two objects
 3. class for sphere?
 4. determination of polyhedrons with vertices => probably need convex hull construction
*/

/**	\brief An abstract struct to define particles shapes 
 *	Shape struct provides a common interfaces for various types of particle shapes, 
 *	such as 2D polygons, 3D polyhedrons, and some special shapes. */

struct Shape {
protected:
	DimensionType Dimension;	/**< space dimension */
	mutable double R_circum,	/**< radius of a sphere circumscribing the shape*/
		R_inscr,				/**< radius of a sphere inscribed inside the shape.*/
		volume,					/**< the volume (in 3D) or area (in 2D) of a particle  */
		area;					/**< the surface area (in 3D) or boundary length (in 2D) of a particle*/
	std::string name;			/**< the name of derived class.*/
	GeometryVector Centroid;	/**< The centroid assuming uniform density. */

	/**	A pure virtual member function to measure the volume of a particle and overwrite this->volume. 	
	 *	Please specifically redefine this function for each particle shape.	 */
	virtual void UpdateVolume() = 0;

	/** A pure virtual member function to measure and overwrite this->area.	 */
	virtual void UpdateSurfaceArea() = 0;

	/** A pure virtual member function to measure and overwrite this->R_circum and this->R_inscr.	 */
	virtual void UpdateR() = 0;

public:
	/** A default constructor	 */
	Shape() {
		this->Dimension = 1;
		this->R_circum = -1;	this->R_inscr = -1.0; this->volume = -1;		this->area = -1;
	};
	
	/** A copy constructor.	 
	 * It only copies dimensions, r_circum, volume, and area. */
	Shape(const Shape & source) {
		this->Dimension = source.Dimension;
		this->R_circum = source.R_circum;	this->volume = source.volume;		this->area = source.area;
		this->R_inscr = source.R_inscr;
	}

	/* Basic interfaces */
	
	/** A member function to return space dimension.
	 * \return space dimension of this particle.	 */
	inline size_t GetDimension() const { return this->Dimension; };
	//virtual void GetVertices(std::vector<GeometryVector> & vertices) const = 0;

	/** A member function to print information of this particle on a screen. */
	void Print() const { Print(std::cout); };
	
	/** A pure virtual function to print information of this particle through ofile.
	 * \param ofile	 output stream. 	 */
	virtual void Print(std::ostream & ofile) const = 0;

	/** A pure virtual member function to move particles by x.  */
	virtual void ShiftBy(const GeometryVector & x) = 0;

	/** A member function to return the volume of this particle.
	 * It asks to compute the volume of this particle, only when it hasn't been computed yet.
	 * \see UpdateVolume()
	 * \return particle volume (in 3D) or the area (in 2D).	 */
	inline double GetVolume() {
		if (this->volume < 0)
			this->UpdateVolume();

		return this->volume;
	}

	/** A member function to return the volume of this particle.
	 * \see GetVolume(); */
	inline double GetVolume() const { return this->volume; };

	/** A member function to return the surface area of this particle.
	 * It asks to compute the surface area of this particle, only when it hasn't been computed yet.
	 * \see UpdateSurfaceArea()
	 * \return the surface area (in 3D) or the boundary length (in 2D).	 */
	inline double GetSurfaceArea() {
		if (this->area < 0)
			this->UpdateSurfaceArea();

		return this->area;
	};

	/** A member function to return the surface area of this particle.
	 * \see GetSurfaceArea()	 */
	inline double GetSurfaceArea() const { return this->area; };

	/** A member function to return the radius of the sphere that circumscribes this particle.
	 * It asks to compute the surface area of this particle, only when it hasn't been computed yet.
	 * \see UpdateR();
	 * \return the maximal distance from the centroid.	 */
	inline double GetMaxRadius() {
		if (this->R_circum < 0) {
			this->UpdateR();
		}
		return this->R_circum;
	};

	/** A member function to return the radius of the sphere that inscribes this particle.
	 * It asks to compute the surface area of this particle, only when it hasn't been computed yet.
	 * \see UpdateR();
	 * \return the least distance from the centroid to its boundary.	 */
	inline double GetInscrRadius() {
		if (this->R_inscr < 0) {
			this->UpdateR();
		}
		return this->R_inscr;
	};

	/** A member function to return the radius of the sphere that circumscribes this particle.
	 * \return the radius = 0 if it is not updated.
	 * \see GetMaxRadius()	 */
	inline double GetMaxRadius() const {
		if (this->R_circum > 0)
			return this->R_circum;
		else {
			return 0.0;
		}
	};

	/** A member function to return the radius of the sphere that inscribes this particle. (CLAUDE)
	 * \return the radius = 0 if it is not updated.
	 * \see GetInscrRadius()	 */
	inline double GetInscrRadius() const {
		if (this->R_inscr > 0)
			return this->R_inscr;
		else {
			return 0.0;
		}
	};
	
	/** A pure virtual member function to compute the center of mass (centroid) of a Shape object. 
	 * \return The centroid.	*/
	virtual GeometryVector GetCenterOfMass() = 0;


	/**	A pure virtual member function to check whether a given point is inside this shape.
	 * \param pos_wrt_centroid	A GeometryVector of a test point in Cartesian coordinates with respect to the particle centroid.
	 * \return True = pos is inside this Shape.	 */
	virtual bool IsPtInside(const GeometryVector & pos_wrt_centroid) { 
		std::cout << "Shape::IsPtInside is not implemented. Please write an implementation in each class." << std::endl; 
		return false;
	}

	/**	A pure virtual member function to compute the closest distance from a given point to the boundary of the shape.
	 * \param pos_wrt_centroid	A GeometryVector of a test point in Cartesian coordinates with respect to the particle centroid.
	 * \param[out] list_dist	A list of distances to sides. list_dist[i] = {ith shortest distance, parameter}.
	 * \param[out] x0	The point on the particle closest to pos_wrt_centroid. 
 	 * \param[out] n_outward	The unit outward normal vector at x0. 
	 * \return distance.	 */
	virtual double GetDistance2Boundary(const GeometryVector & pos_wrt_centroid, std::vector<GeometryVector> & list_dist, GeometryVector & x0, GeometryVector & n_outward ) const { 
		std::cout << "Shape::GetDistance2Boundary is not implemented. Please write an implementation in each class." << std::endl; 
		return 0.0;
	}
	
	/**	A pure virtual member function to compute the closest distance from a given point to the boundary of the shape.
	 * \param pos_wrt_centroid	A GeometryVector of a test point in Cartesian coordinates with respect to the particle centroid.
	 * \return distance.	 */
	virtual double GetDistance2Boundary(const GeometryVector & pos_wrt_centroid) const { 
		std::cout << "Shape::GetDistance2Boundary is not implemented. Please write an implementation in each class." << std::endl; 
		return 0.0;
	}

	/** A pure virtual member function to copy this Shape class.
	 * It will be useful when the object is stored in a pointer. 
	 * \return A shared pointer of this Shape object, which includes identical information.
				Since this is a "copy," we don't need to worry about security. */
	virtual std::shared_ptr<Shape> clone() = 0;

	//virtual std::shared_ptr<Shape> copy() const = 0;
	/*Features in future
	//TODO
	virtual bool Overlap(const GeometryVector & adjParticle, const GeometryVector & Coordinates) = 0;
	virtual bool VoxelOverlap( ...)
	virtual void ProcessSplittedVoxels(...)
	*/


/* Transformations 
	Rotate	rotate vertices around its center
	Scale	expand/shrink particle shape by "factor"
*/
	/**	A pure virtual member function to rotate this particle.
	 * It rotates Cartesian coordinates to from the basis (1,0,0), (0,1,0), (0,0,1) to a new basis  
	 * \param newbasis	Orthonormal basis vectors of a new Cartesian Coordinates.	 */
	virtual void Rotate(const std::vector<GeometryVector> & newbasis) = 0;

	/** A member function to rotate this particle into a random orientation.
	 * \param rng A random generator to make random orthonormal basis vectors
	 * \see Rotate(const std::vector<GeometryVector> & newbasis);	 */
	void Rotate(RandomGenerator & rng);
	
	/** A member function to resize this particle to a prescribe volume.
	 * \param newVolume	Desired volume of this particle.	 */
	void Scale2Vol(double newVolume);

	/** A member function to resize this particle at a linear scale.
	 * \param factor	The linear scale factor at which this particle is resized.
	 * \see Scale2Vol(double newVolume)	 */
	virtual void Scale(double factor);
	
	/** A pure virtual member function to compute a Fourier transform of this particle.
	 * Fourier transform is computed by the following convention: \int \d{\vect{r}} e^{i k r} f(r).
	 * Form factor is computed with respect to the centroid of this particle.
	 * \param k	A wavevector at which the form factor (Fourier transform) is computed.
	 * \return Complex Fourier transform of this particle	 */
	virtual std::complex<double> GetFormFactor(const GeometryVector & k) = 0;
};
/** 
 * A function to check the orthonormality of a given basis vectors.
 * \param dimension	Space dimension
 * \param newbasis	d GeometryVectors in d dimensions
 * \return true = basis vectors are orthonormal.
 */
bool CheckNewCoordinates(int dimension, const std::vector<GeometryVector> & newbasis);

/** A shift operator << for the Shape classes. 
 * \param os	An output stream.
 * \param shape	A Shape object. */
inline std::ostream & operator << (std::ostream & os, const Shape & shape) {
	shape.Print(os);
	return os;
}

/** \brief A child struct from the Shape struct to define a planar figure in 3D.
 * This is a polygonal plane in 3D Euclidean space and thus have orientations.
 * This is an auxiliary struct to define a 2D polygon and facets of a 3D polyhedron. However, we seldom use this struct alone. 
 * \see Shape
 */
struct PlanarFigure : public Shape {
	friend struct Polygon;
	friend struct Polyhedron;
protected:
	std::vector<GeometryVector *> V;	///< Pointers of Cartesian coordinates of vertices.
										///< Vertex indices that consist of this planar figure
										///< vert_idx.size() <=  num_prts
										///< V[vert_idx[i]] becomes a cyclic chain of vertices.
	size_t num_Vs;						///< the number of vertices
	double r_perp;						///< distance from the origin to the plane
	mutable GeometryVector normal;		///< normal vector
	mutable bool Trust_V;				///< false = check whether vertices on a plane.
	mutable bool update_ER;				///< false = need to update E and R vectors
										
	std::vector<GeometryVector> E;		///< relative positions of a vertex with respec to its neighboring vertex 
										///< auxiliary vectors to compute the form factor
	std::vector<GeometryVector> R;		///< the mean of two neighboring vertices. 
										///< auxiliary vectors to compute the form factor
	
	/** 
	 * A member function to overwite this->volume with the area of this facet.
	 * \see Shape::UpdateVolume()
	 */
	void UpdateVolume();

	/** 
	 * A member function to overwrite this->area with the total side length of this facet.
	 * \see Shape::UpdateSurfaceArea();
	 */
	void UpdateSurfaceArea();

	/** A member function to overwrite this->R_circum with the maximal distance from the origin to vertcies.
	 * It doesn't udpate this->R_inscr.
	 * \see Shape::UpdateR()
	 */
	void UpdateR();	

	/**
	 * A member function to check whether all vertices are on the same plane, and in a cyclic order.
	 * It also updates the normal vector.
	 */
	void CheckVertices();		

	/** 
	 *	A member function to update E_i and R_i vectors
	 */
	void Update_ER_vectors();


public:
	/* Constructors */
	//PlanarFigure(const PlanarFigure & source);					//copy
	/** A default constructor */
	PlanarFigure() : Shape() { this->r_perp = 0.0; };
	
	/** A constructor when all points are on this facet.
	 *  \param V	pointers of GeometryVectors of vertices in d dimensions
	 *	\param num_vertices	The number of vertices.
  	 */
	PlanarFigure(GeometryVector * V, size_t num_vertices);			
	
	/** A constructor when a portion of points forms this facet 
	 * \param V	pointers of GeometryVectors of vertices in d dimension.
	 * \param vert_idx	Index of vertices that form this facet. V[vert_idx[0]], ... V[vert_idx[n-1]] should form a cyclic chain of vertices on a plane.
	 */
	PlanarFigure(GeometryVector * V, std::vector<size_t> vert_idx);	

	/* Basic Interfaces */
	void Print(std::ostream & ofile) const;	///< \see Shape::Print(std::ostream & ofile)

	/**
	 *	A member function to obtain Cartesian coordinates of vertices that form this facet.
	 * \param[out]	Verts	Verts[i] = a GeometryVector of ith vertex of this facet.
	 */
	void GetVertices(std::vector<GeometryVector> & Verts) const;
	
	/** 
	 * A member function to compute the centor of mass of this facet with respect to the origin
	 * \return the center of mass of this facet.
	 */
	GeometryVector GetCenterOfMass();

	/** 
	 *	A member function to rotate this facet (it olny rotates the normal vector).
	 * It only works in 3D.
	 */
	void Rotate(const std::vector<GeometryVector> & newbasis);

	/** A member function to resize this facet at a linear scale
	 * This DO NOT rescale vertices.
	 * \param factor	the linear factor
	 */
	void Scale(double factor);		

	/** A member function to compute the Form factor of this facet.
	 * \param k	A wavevector at which the form factor is computed
	 * \return A complex Fourier component at k.
	*/
	std::complex<double> GetFormFactor(const GeometryVector & k);

	/** Please don't use this member funciton. */
	virtual std::shared_ptr<Shape> clone() { return nullptr; }

	/** 
	 * A memebr function to move the position of a PlanarFigure object.
	 * No implementation!!
	*/
	void ShiftBy(const GeometryVector & x){}
	// Inherited
private:

	/** 
	 * A member function to compute the form factor of this facet.
	 * It uses an analytic expression, which can badly behave at small wavevectors.
	 * \see GetFormFactor(const GeometryVector & k)
	 * \param k	A wavevector at which the form factor is computed
	 * \return A complex Fourier component at k.
	 */
	std::complex<double> FormFactor_analytic(const GeometryVector &k);

	/** 
	 * A member function to compute the form factor of this facet.
	 * It uses a series expansion formula, which can well behave at small wavevectors.
	 * \see GetFormFactor(const GeometryVector & k)
	 * \param k	A wavevector at which the form factor is computed
	 * \return A complex Fourier component at k.
	 */
	std::complex<double> FormFactor_series(const GeometryVector &k);
};


/** \brief A struct to define 2D general polygon.
 * It consists of vertices and an object of PlanarFigure 
 * The origin is set to be the center of mass
 */
struct Polygon : public Shape {
protected:
	
	std::vector<GeometryVector> Vertices;	///< Cartesian coordinates of vertices with respect to Centroid.
	PlanarFigure Gamma;						///< A facet.

	/** A member function to overwrite R_circum with the maximal distance from the centroid to vertices. */
	void UpdateR();					

	/** A member function to overwrite this->volume with the area of this polygon. */
	void UpdateVolume();

	/** A member function to overwrite this->area with the total side length of this polygon.	 */
	void UpdateSurfaceArea();

	/** A common member function to initialize the object. */
	void init(const GeometryVector * vertices, size_t num_vertices);


public:
	Polygon() :Shape::Shape() { this->Dimension = 2; name = "Polygon"; }; ///< A default constructor

	/** A copy constructor.  */
	Polygon(const Polygon & source);

	/** A constructor by vertices that form this facet.
	 * \see Polygon(const GeometryVector *, size_t);
	 * \param vertices	GeometryVectors in the cyclic order	 */
	Polygon(const std::vector<GeometryVector> & vertices);
	
	/** A constructor by vertices that form this facet.
	 * Vertices are automatically translated in such a way that the origin coincides with the center of the polygon.
	 * \param vertices	An array of GeometryVectors that will form this facet. They should be a cyclic order
	 * \param num_vertices	The number of vertices 	 */
	Polygon(const GeometryVector * vertices, size_t num_vertices);
	
	/** A constructor by string input stream. 
	 * \param ifile string input stream.*/
	Polygon(std::istream & ifile);


	/** A member function to print the number of vertices and Cartessian coordinates of vertices.
	 * \see Shape::Print(std::ostream &)	 */
	void Print(std::ostream & ofile) const {
		ofile << "[" << name << "]\n";
		std::vector<GeometryVector> vs;
		this->GetVertices(vs);
		ofile << "#Vertices\t"<<vs.size()<<"\n";
		for (int i=0; i<vs.size(); i++)
			ofile << vs[i];

	}

	/** A member function to make the indices of vertices "positive." */
	void SetPositiveOrientation();

	/// A member function = Print(std::cout) 
	/// \see Polygon::Print(std::ostream)
	void Print() const { this->Print(std::cout); };

	/** A member function to get all vertices of this polygon
	 * \param[out]	vertices	GeometryVectors of all vertices.	 */
	void GetVertices(std::vector<GeometryVector> & vertices) const;

	/**	A member function to compute the center of mass of this polygon
	 * \return the centroid.	 */
	GeometryVector GetCenterOfMass();

	/**	A virtual member function to check whether a given point is inside this Polygon.
	 * \param pos_wrt_centroid	A GeometryVector of a test point in Cartesian coordinates with respect to the particle centroid.
	 * \return True = pos is inside this Shape.	 */
	virtual bool IsPtInside(const GeometryVector & pos_wrt_centroid);

	/**	A member function to compute the closest distance from a given point to the side of this polygon.
	 * \param pos_wrt_centroid	A GeometryVector of a test point in Cartesian coordinates with respect to the particle centroid.
	 * \return distance.	 */
	virtual double GetDistance2Boundary(const GeometryVector & pos_wrt_centroid) const;

	/**	A member function to compute the closest distance from a given point to the side of this polygon.
	 * \param pos_wrt_centroid	A GeometryVector of a test point in Cartesian coordinates with respect to the particle centroid.
	 * \param[out] list_dist	A list of distances in the ascending order.
 	 * \param[out] x0	The point on the particle closest to pos_wrt_centroid. 
 	 * \param[out] n_outward	The unit outward normal vector at x0. 
	 * \return distance.	 */
	virtual double GetDistance2Boundary(const GeometryVector & pos_wrt_centroid, std::vector<GeometryVector> & list_dist, GeometryVector & x0, GeometryVector & n_outward) const;

	/** \return A shared pointer of an identical copy of this Polygon object. */
	virtual std::shared_ptr<Shape> clone() {
		return std::shared_ptr<Shape>(new Polygon(*this));
	}

	/** A member function to move all vertices by x.
	 * E and R vectors in this->Gamma are also updated.
	 * \param x	A GeometryVector by which all vertices are translated.	 */
	void ShiftBy(const GeometryVector & x);

	/// \see Shape::Rotate(RandomGenerator)
	void Rotate(RandomGenerator & rng) { Shape::Rotate(rng); };
	/// \see Shape::Rotate(RandomGenerator)
	void Rotate(const std::vector<GeometryVector> & newbasis);
	/// \see Shape::Scale(double)
	virtual void Scale(double factor);

	/** A member function to compute the form factor of this 2D polygon with respect to its centroid at a wavevector k.
	 * Currently, it directly recalls the form factor of an object Gamma. 
	 * \param k	A wavevector at which we compute the form factor
	 * \return a complex form factor.	 */
	inline std::complex<double> GetFormFactor(const GeometryVector & k) { return this->Gamma.GetFormFactor(k); }
};

/** Read a Polygon from a string input stream. (CLAUDE)
 * \param ifile	A string input stream holding the number of vertices and their coordinates.
 * \return The Polygon reconstructed from the stream. */
Polygon ReadPolygon(std::istream & ifile);

/** \brief A struct to define 3D general polyhedron.
 * It consists of vertices, lists of vertices that form each facet, and objects of PlanarFigure, each of which describes the facet.
 * The origin is set to be the center of mass. */
struct Polyhedron : public Shape{
protected:
	std::vector<GeometryVector> Vertices;			///< Cartesian coordinates of all vertices that form this Polyhedron.
	std::vector<PlanarFigure> Gamma;				///< Gamma[i] describes the ith facet.
	std::vector<std::vector<size_t>> Vertex_idx;	///< Vertex_idx[i][j]: the index of the jth vertex in Gamma[i]
	

	/** A member function to overwrite this->R_circum = the largest distance from the centroid to vertices,
	 *								and this->R_inscr = the smallest r_perp.	 */
	void UpdateR();	
	
	/**	A member function to overwrite this->volume with the volume of this polyehdron.	 */
	void UpdateVolume();

	/** A member function to overwrite this->area with the surface area of this polyhedron.	 */
	void UpdateSurfaceArea();

	/**	A member function to check whether this object satisfies the Euler identity.
	 * \return true = satisfied.	 */
	bool CheckEulerIdentity() const;

	/** A member function to check whether this object is enclosed.
	 *	It adds the normal vectors of all facets, which is supposed to be zero for an enclosed polyhedron.
	 * \return true = Summation of all normal vectors = 0.	 */
	bool CheckClosedSuface();

	/** A member function to check whether all normal vectors outwards.
	 *	If not, it changes directions of normal vectors.	 */
	void CheckOrientations();						


public:

	/**	A default constructor for Polyhedron. */
	Polyhedron() : Shape::Shape() { this->Dimension = 3; name = "Polyhedron"; };

	/** A copy constructor for a Polyhedron.
	 *  It blidnly copies the source.
	 * \param source	A Polyhedron obejct. */
	Polyhedron(const Polyhedron & source) : Polyhedron::Polyhedron(source.Vertices, source.Vertex_idx) {
		this->volume = source.volume;	this->area = source.area;	this->R_circum = source.R_circum; this->R_inscr = source.R_inscr;
		//this->Centroid = GeometryVector((int)0);
	}

	/** A constructor of a Polyhedron object by vertices and vertex indices that form each facet.
	 *	This constructor checks that this polyhedron can be exactly constructed.
	 * \param vertices	Cartesian coordinates of all vertices.
	 * \param idx	vertices[idx[i][j]] = the jth vertex of the ith facet of this Polyhedron. */
	Polyhedron(const std::vector<GeometryVector> & vertices, 
		const std::vector<std::vector<size_t>> & idx);

	/** A constructor by a string input stream. 
	 * \param ifile string input stream.*/
	Polyhedron(std::istream & ifile);

	/* Can we determine faces only with vertices?? */
	
	/** A member function to print the number of vertices, cartesian coordinates of all vertices, the number of faces, and vertex indices that form each faces.
	 *  See example: a cube
	 *  [Polyhedron]
	 *	# Vertices:	8
	 *	-1,	-1,	-1,	0
	 *	-1,	-1,  1,	0
	 *	....
	 *	# Faces: 6
	 *  0	1	2	3
	 *	0	1	...

	 *	\param ofile	output stream
	 *	\see Shape::Print(std::ostream &)
	 */
	void Print(std::ostream & ofile) const;

	/** A member function to print on a screen the number of vertices, cartesian coordinates of all vertices, the number of faces, and vertex indices that form each faces.
	 *	Polyhedron::Print(std::cout)	 */
	void Print() const { this->Print(std::cout); };

	/** A member function to obtain cartesian coordinates fo all vertices, and vertex indices that can reconstruct this Polyhedron.
	 *	\param[out]	vertices	Cartesian coordinates of all vertices.
	 *	\param[out]	idx[i][j]	The jth vertex index that forms the ith face.	 */
	void GetVertices(std::vector<GeometryVector> & vertices, 
		std::vector<std::vector<size_t>> & idx) const;

	/** A member function to compute the center of mass.
	 *	\return Cartesian coordinates of the centroid with respect to the origin. */
	GeometryVector GetCenterOfMass();

	/** A member function to randomly rotate this Polyhedron
	 * \param rng	A random generator to construct a new set of basis vectors.	 */
	void Rotate(RandomGenerator & rng) { Shape::Rotate(rng); };

	/** A member function to rotate this Polyhedron to a new cartesian coordiantes.
	 *	\param newbasis	A new set of basis vectors in cartesian coordinates.	 */
	void Rotate(const std::vector<GeometryVector> & newbasis);

	/** A member function to resize this Polyhedron at a linear scale "factor".
	 * \param factor	A linear factor at which this Polyhedron is rescaled.	 */
	virtual void Scale(double factor);

	/** A member function to translate all cartesian coordinates of vertices and associated E and R vectors by x.	 
	 * \param x	A GeometryVector by which this Polyhedron is translated. */
	void ShiftBy(const GeometryVector & x);

	/**	A virtual member function to check whether a given point is inside this Polygon.
	 * \param pos_wrt_centroid	A GeometryVector of a test point in Cartesian coordinates with respect to the particle centroid.
	 * \return True = pos is inside this Shape.	 */
	virtual bool IsPtInside(const GeometryVector & pos_wrt_centroid);

	/**	A pure virtual member function to compute the closest distance from a given point to the facets of this polygon.
	 * \param pos_wrt_centroid	A GeometryVector of a test point in Cartesian coordinates with respect to the particle centroid.
	 * \return distance.	 */
	virtual double GetDistance2Boundary(const GeometryVector & pos_wrt_centroid) const;

	/**	A pure virtual member function to compute the closest distance from a given point to the facets of this polygon.
	 * \param pos_wrt_centroid	A GeometryVector of a test point in Cartesian coordinates with respect to the particle centroid.
	 * \param[out] list_dist	A list of distances to the point in the ascending order.
	 * \param[out] x0	The point on the particle closest to pos_wrt_centroid. 
 	 * \param[out] n_outward	The unit outward normal vector at x0. 
	 * \return distance.	 */
	virtual double GetDistance2Boundary(const GeometryVector & pos_wrt_centroid, std::vector<GeometryVector> & list_dist, GeometryVector & x0, GeometryVector & n_outward) const;

	/** \return A shared pointer of an identical copy of this Polyhedron object. (CLAUDE) */
	virtual std::shared_ptr<Shape> clone() {
		return std::shared_ptr<Shape>( new Polyhedron (*this));
	}

	/** A member function to compute the form factor of this 3D polyhedron at a wavevector k.
  	 * Currently, it directly recalls the form factor.
	 * \param k	A wavevector at which we compute the form factor
	 * \return a complex form factor	 */
	std::complex<double> GetFormFactor(const GeometryVector & k);
private:

	/** A member function to compute the form factor of this 3D polyherdon at a wavevector k.
	 *	It uses an analytic formula, which may give incorrect results for small wavevectors due to the singularity at the origin.
	 * \param k	A wavevector at which we compute the form factor
	 * \return a complex form factor	 */
	std::complex<double> FormFactor_analytic(const GeometryVector &k);

	/**	A member function to compute the form factor of this 3D polyherdon at a wavevector k.
	 *	It uses a series expansion of the form factor which becomes expensive at large wavevectors.
	 *	It is not implemented yet.
	 * \param k	A wavevector at which we compute the form factor
	 * \return a complex form factor	 */
	std::complex<double> FormFactor_series(const GeometryVector &k);
};

//TODO:
//implement pre-defined polygons and polyhedrons, e.g., square, cube...

///* Parallelepiped shape*/
//struct Parallelepiped : public Shape {
//protected:
//	std::vector<GeometryVector> a; /* basis vector for sides*/
//	std::vector<GeometryVector> b; /* reciprocal vectors */
//
//public:
//	// vertices = 0.5*(\pm sides[0] + \pm sides[1] + ...)
//	Parallelepiped(const Parallelepiped & source);
//	Parallelepiped(DimensionType dimension, const std::vector<GeometryVector> & sides);
//	Parallelepiped(DimensionType dimension, double L);	/* L = side length*/
//	void GetVertices(std::vector<GeometryVector> & vertices) const;
//	void Print() const;	/* Print axes */
//	virtual double GetVolume() const;
//	virtual double GetSurfaceArea() const;
//
//	virtual void Rotate(const std::vector<GeometryVector> & newbasis);
//	void Rotate(RandomGenerator & rng) { Shape::Rotate(rng); };
//	virtual void Scale(double factor);
//	virtual void UpdateR() const;
//
//	virtual const GeometryVector & GetBasis(size_t i) const { return this->a[i]; }
////	void Rotate(std::vector<GeometryVector> & newcoordinatebasis);
////	void Rotate(RandomGenerator & rng);
////	void 
//};
//
///* Cubic */
//struct Cube : public Parallelepiped{
//protected:
//	double L, halfL;
//public:
//	/* Generate a cube with sides (L,0,0), (0,L,0), (0,0,L) in d=3*/
//	Cube(const Cube & source);
//	Cube(DimensionType dimension, double L);
//
//	double GetSideLength() const { return this->L; };
//	double GetHalfLength() const { return this->halfL; };
//};
//


#endif