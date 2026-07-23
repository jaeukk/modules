/**
*	Author	: Jaeuk Kim
*	Email	: phy000.kim@gmail.com
*	Date	: April. 2019 */

/** \file PeriodicNetworks.cpp
* \brief Function implementations for "PeriodicNetworks" and a CMD interface. */



#include	"PeriodicNetworks.h"
#include	<StructureFactor.h>
#include	<GenerateConfigs.h>
#include	<NumberVariance.h>
#include	<SpectralDensity_Dispersion.h>
#include	<SpectralDensity.h>
#include	<VolumeFractionVariance.h>



void GetPeriodicPackings(std::istream & ifile, std::ostream & ofile, SpherePacking &sphC, std::string & filename) {
	double	rho = 0,	/**< number density */
		D = 0,			/**< nearest-neighbor distances*/
		a = 0;			/**< lattice constant */

	std::string NetworkName, Lengthscale;

	ofile << "Network name = ";	ifile >> NetworkName;
	ofile << "Lengthscale (D/ rho/ a) = "; ifile >> Lengthscale;

	if (Lengthscale.compare("D") == 0) {
		ofile << "Nearest Neighbor distance, D = ";	ifile >> D;
	}
	else if (Lengthscale.compare("rho") == 0) {
		ofile << "Number density, rho = ";	ifile >> rho;
	}
	else if (Lengthscale.compare("a") == 0) {
		ofile << "Lattice constant = ";	ifile >> a;
	}
	else {
		ofile << "Choose a correct length scale !\n";
		return;
	}

	LatticeGen lattice;
	double factor = 1.0;

	/*Define periodic network*/
	if (NetworkName.compare("custom") == 0) {
		ofile << "custom: is not implemented yet.\n";
		//TODO: write a code for generating a custom network.
	}
	else {
		lattice = LatticeGen(NetworkName);
		lattice.GenerateP(sphC);

		double NumDensity = sphC.NumParticle() / sphC.PeriodicVolume();
		if (Lengthscale.compare("D") == 0) {
			factor = D / sphC.NearestParticleDistance(0);
		}
		else if (Lengthscale.compare("rho") == 0) {
			factor = pow(rho / NumDensity, 1.0 / sphC.GetDimension());
		}
		else if (Lengthscale.compare("a") == 0) {
			factor = 1.0 / sqrt(sphC.GetBasisVector(0).Modulus2());
		}
		sphC.Rescale(factor);
	}
	filename = NetworkName + "_spheres";
}

void GetPeriodicPackings(const std::string & NetworkName, const std::string & Lengthscale, double val, SpherePacking & sphC, std::string & filename) {
	double	rho = 0,	/**< number density */
		D = 0,			/**< nearest-neighbor distances*/
		a = 0;			/**< lattice constant */

	//std::string NetworkName, Lengthscale;

	std::cout << "Network name = " << NetworkName <<"\n";
	std::cout << "Lengthscale (D/ rho/ a) = " << Lengthscale <<"\n";

	if (Lengthscale.compare("D") == 0) {
		std::cout << "Nearest Neighbor distance, D = ";  
		D = val;
	}
	else if (Lengthscale.compare("rho") == 0) {
		std::cout << "Number density, rho = ";	
		rho = val;
	}
	else if (Lengthscale.compare("a") == 0) {
		std::cout << "Lattice constant = ";	
		a = val;
	}
	else {
		std::cout << "Choose a correct length scale !\n";
		return;
	}

	LatticeGen lattice;
	double factor = 1.0;

	/*Define periodic network*/
	if (NetworkName.compare("custom") == 0) {
		std::cout << "custom: is not implemented yet.\n";
		//TODO: write a code for generating a custom network.
	}
	else {
		lattice = LatticeGen(NetworkName);
		lattice.GenerateP(sphC);

		double NumDensity = sphC.NumParticle() / sphC.PeriodicVolume();
		if (Lengthscale.compare("D") == 0) {
			factor = D / sphC.NearestParticleDistance(0);
		}
		else if (Lengthscale.compare("rho") == 0) {
			factor = pow(rho / NumDensity, 1.0 / sphC.GetDimension());
		}
		else if (Lengthscale.compare("a") == 0) {
			factor = a / sqrt(sphC.GetBasisVector(0).Modulus2());
		}
		sphC.Rescale(factor);
	}
	filename = NetworkName + "_spheres";
}

void GetPeriodicNetworks(const std::string & NetworkName, const std::string & Lengthscale, double val, Dispersion & polyC, std::string & filename)
{
	double	rho = 0,	/**< number density */
		D = 0,			/**< nearest-neighbor distances*/
		a = 0;			/**< lattice constant */

	

	std::cout << "Network name = " << NetworkName;
	std::cout << "Lengthscale (D/ rho/ a) = " << Lengthscale;

	if (Lengthscale.compare("D") == 0) {
		D = val;
		std::cout << "Nearest Neighbor distance, D = " << D;	
	}
	else if (Lengthscale.compare("rho") == 0) {
		rho = val;
		std::cout << "Number density, rho = " << rho;
	}
	else if (Lengthscale.compare("a") == 0) {
		a = val;
		std::cout << "Lattice constant = " << a;
	}
	else {
		std::cout << "Choose a correct length scale !\n";
		return;
	}

	LatticeGen lattice;
	Configuration c;
	double factor = 1.0;

	/*Define periodic network*/
	if (NetworkName.compare("custom") == 0) {
		std::cout << "custom: is not implemented yet.\n";
		//TODO: input for custom periodic networks.
	}
	else {
		/*=========2D========*/
		if (NetworkName.compare("SquareNetwork") == 0) {
			lattice = LatticeGen("Square");
			lattice.GenerateC(c);
			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1);
			std::vector<GeometryVector> vertices;
			vertices.emplace_back(0.5*(v1 + v2));
			vertices.emplace_back(0.5*((-1.0)*v1 + v2));
			vertices.emplace_back(0.5*((-1.0)*v1 + (-1.0)*v2));
			vertices.emplace_back(0.5*(v1 - v2));

			polyC = Dispersion(c, Polygon(vertices).clone());
		}
		else if (NetworkName.compare("OctagonalNetwork") == 0) {
			lattice = LatticeGen("Square");
			lattice.GenerateC(c);
			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1);
			std::vector<GeometryVector> vertices;
			double f = 1.0 / (2.0 + sqrt(2.0));
			/*The first particle (square)*/
			vertices.emplace_back(f*v1);
			vertices.emplace_back(f*v2);
			vertices.emplace_back(-f*v1);
			vertices.emplace_back(-f*v2);
			polyC = Dispersion(c, Polygon(vertices).clone());
			polyC.MoveParticle(0, GeometryVector(0.0, 0.0));

			/*The second particle (octagon)*/
			GeometryVector x0 = 0.5*(v1 + v2);
			vertices[0] = f*v1 - x0;
			vertices[1] = (1.0-f)*v1 - x0;
			vertices[2] = v1 + f*v2 - x0;
			vertices[3] = v1 + (1.0-f)*v2 - x0;
			vertices.emplace_back((1.0-f)*v1 + v2 - x0);
			vertices.emplace_back(f*v1 + v2 - x0);
			vertices.emplace_back((1.0-f)*v2 - x0);
			vertices.emplace_back(f*v2 - x0);
			polyC.Insert(Polygon(vertices).clone(), GeometryVector(0.5, 0.5));

		}
		else if (NetworkName.compare("RhombicNetwork") == 0) {
			lattice = LatticeGen("Triangular");
			lattice.GenerateC(c);
			/*The first particle*/
			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1),
				x0 = 0.5*(v1 + v2);
			std::vector<GeometryVector> vertices;
			vertices.emplace_back(0.0*v1 - x0);
			vertices.emplace_back(v1 - x0);
			vertices.emplace_back(v1 + v2 - x0);
			vertices.emplace_back(v2 - x0);
			polyC = Dispersion(c, Polygon(vertices).clone());
			polyC.MoveParticle(0, GeometryVector(0.5, 0.5));
		}
		else if (NetworkName.compare("TriangularNetwork") == 0) {
			lattice = LatticeGen("Triangular");
			lattice.GenerateC(c);
			/*The first particle*/
			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1);
			std::vector<GeometryVector> vertices;
			vertices.emplace_back(v1);
			vertices.emplace_back(v2);
			vertices.emplace_back(0.0*v1);
			polyC = Dispersion(c, Polygon(vertices).clone());
			polyC.MoveParticle(0, GeometryVector(1.0 / 3.0, 1.0 / 3.0));

			/*The second particle*/
			vertices[0] = v1;
			vertices[1] = v1 + v2;
			vertices[2] = v2;
			polyC.Insert(Polygon(vertices).clone(), GeometryVector(2.0 / 3.0, 2.0 / 3.0));
		}
		else if (NetworkName.compare("HoneycombNetwork") == 0) {
			lattice = LatticeGen("Triangular");
			lattice.GenerateC(c);

			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1);
			std::vector<GeometryVector> vertices;
			vertices.emplace_back(1.0 / 3.0*(v1 + v2));
			vertices.emplace_back(1.0 / 3.0*((-1)*v1 + 2 * v2));
			vertices.emplace_back(1.0 / 3.0*((-2)*v1 + 1 * v2));
			vertices.emplace_back(1.0 / 3.0*((-1)*v1 - 1 * v2));
			vertices.emplace_back(1.0 / 3.0*((1)*v1 - 2 * v2));
			vertices.emplace_back(1.0 / 3.0*((2)*v1 - 1 * v2));

			polyC = Dispersion(c, Polygon(vertices).clone());
		}
		else if (NetworkName.compare("KagomeNetwork") == 0) {
			lattice = LatticeGen("Triangular");
			lattice.GenerateC(c);
			/*The first particle (triangle)*/
			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1);
			std::vector<GeometryVector> vertices;
			GeometryVector x0 = 1.0 / 6.0*(v1 + v2); //center of mass...
			vertices.emplace_back(0.5*v1 - x0);
			vertices.emplace_back(0.5*v2 - x0);
			vertices.emplace_back(0.0*v1 - x0);
			polyC = Dispersion(c, Polygon(vertices).clone());
			polyC.MoveParticle(0, GeometryVector(1.0 / 6.0, 1.0 / 6.0));

			/*The second particle (triangle) */
			x0 = 5.0 / 6.0*(v1 + v2);
			vertices[0] = v1 + 0.5*v2 - x0;
			vertices[1] = v1 + v2 - x0;
			vertices[2] = 0.5*v1 + v2 - x0;
			polyC.Insert(Polygon(vertices).clone(), GeometryVector(5.0 / 6.0, 5.0 / 6.0));

			/*The third particle (hexagon) */
			x0 = 0.5*(v1 + v2);
			vertices[0] = 0.5*v1 - x0;
			vertices[1] = v1 - x0;
			vertices[2] = v1 + 0.5*v2 - x0;
			vertices.emplace_back(v2 + 0.5*v1 - x0);
			vertices.emplace_back(v2 - x0);
			vertices.emplace_back(0.5*v2 - x0);
			polyC.Insert(Polygon(vertices).clone(), GeometryVector(0.5, 0.5));
		}
		/* ==========3D======== */
		else if (NetworkName.compare("CubicNetwork") == 0) {
			lattice = LatticeGen("Cubic");
			lattice.GenerateC(c);

			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1), v3 = c.GetBasisVector(2),
				x0 = 0.5*(v1 + v2 + v3);
			std::vector<GeometryVector> vertices;	std::vector<std::vector<size_t>> index;
			/*0*/	vertices.emplace_back(0.0*v1 - x0);
			/*1*/	vertices.emplace_back(v2 - x0);
			/*2*/	vertices.emplace_back(v1 + v2 - x0);
			/*3*/	vertices.emplace_back(v1 - x0);
			/*4*/	vertices.emplace_back(v3 - x0);
			/*5*/	vertices.emplace_back(v1 + v3 - x0);
			/*6*/	vertices.emplace_back(v2 + v3 - x0);
			/*7*/	vertices.emplace_back(v1 + v2 + v3 - x0);
			{
				index.emplace_back(std::vector<size_t> { 0, 1, 2, 3 });
				index.emplace_back(std::vector<size_t> { 0, 3, 5, 4 });
				index.emplace_back(std::vector<size_t> { 0, 4, 6, 1 });
				index.emplace_back(std::vector<size_t> { 1, 6, 7, 2 });
				index.emplace_back(std::vector<size_t> { 7, 5, 3, 2 });
				index.emplace_back(std::vector<size_t> { 4, 5, 7, 6 });
			}

			polyC = Dispersion(c, Polyhedron(vertices, index).clone());
		}
		else if (NetworkName.compare("FCCNetwork") == 0) {
			/* Particle = unit cell */
			lattice = LatticeGen("FCC");
			lattice.GenerateC(c);

			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1), v3 = c.GetBasisVector(2),
				x0 = 0.5*(v1 + v2 + v3);
			std::vector<GeometryVector> vertices;	std::vector<std::vector<size_t>> index;

			/*0*/	vertices.emplace_back(0.0*v1 - x0);
			/*1*/	vertices.emplace_back(v2 - x0);
			/*2*/	vertices.emplace_back(v1 + v2 - x0);
			/*3*/	vertices.emplace_back(v1 - x0);
			/*4*/	vertices.emplace_back(v3 - x0);
			/*5*/	vertices.emplace_back(v1 + v3 - x0);
			/*6*/	vertices.emplace_back(v2 + v3 - x0);
			/*7*/	vertices.emplace_back(v1 + v2 + v3 - x0);
			{
				index.emplace_back(std::vector<size_t> { 0, 1, 2, 3 });
				index.emplace_back(std::vector<size_t> { 0, 3, 5, 4 });
				index.emplace_back(std::vector<size_t> { 0, 4, 6, 1 });
				index.emplace_back(std::vector<size_t> { 1, 6, 7, 2 });
				index.emplace_back(std::vector<size_t> { 7, 5, 3, 2 });
				index.emplace_back(std::vector<size_t> { 4, 5, 7, 6 });
			}

			polyC = Dispersion(c, Polyhedron(vertices, index).clone());
		}
		else if (NetworkName.compare("BCCNetwork") == 0) {
			/* particle = unit cell */
			lattice = LatticeGen("BCC");
			lattice.GenerateC(c);

			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1), v3 = c.GetBasisVector(2),
				x0 = 0.5*(v1 + v2 + v3);
			std::vector<GeometryVector> vertices;	std::vector<std::vector<size_t>> index;

			/*0*/	vertices.emplace_back(0.0*v1 - x0);
			/*1*/	vertices.emplace_back(v2 - x0);
			/*2*/	vertices.emplace_back(v1 + v2 - x0);
			/*3*/	vertices.emplace_back(v1 - x0);
			/*4*/	vertices.emplace_back(v3 - x0);
			/*5*/	vertices.emplace_back(v1 + v3 - x0);
			/*6*/	vertices.emplace_back(v2 + v3 - x0);
			/*7*/	vertices.emplace_back(v1 + v2 + v3 - x0);
			{
				index.emplace_back(std::vector<size_t> { 0, 1, 2, 3 });
				index.emplace_back(std::vector<size_t> { 0, 3, 5, 4 });
				index.emplace_back(std::vector<size_t> { 0, 4, 6, 1 });
				index.emplace_back(std::vector<size_t> { 1, 6, 7, 2 });
				index.emplace_back(std::vector<size_t> { 7, 5, 3, 2 });
				index.emplace_back(std::vector<size_t> { 4, 5, 7, 6 });
			}

			polyC = Dispersion(c, Polyhedron(vertices, index).clone());
		}
		else if (NetworkName.compare("HCPNetwork") == 0) {
			/* particle = unit cell */
			//TODO: complete it
			lattice = LatticeGen("HCP");
			lattice.GenerateC(c);
		}
		else if (NetworkName.compare("FCCNewtork_WS") == 0) {
			/* Particle = Wigner-Seitz cell */
			//TODO: complete
			//lattice = LatticeGen("FCC");
			//lattice.GenerateP(sphC);

			//GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1), v3 = c.GetBasisVector(2);
			//std::vector<GeometryVector> vertices;	std::vector<std::vector<size_t>> index;

		}
		else if (NetworkName.compare("BCCNetwork_WS") == 0) {
			//TODO: complete it
			lattice = LatticeGen("BCC");
			lattice.GenerateC(c);
		}
		else if (NetworkName.compare("HCPNetwork_WS") == 0) {
			//TODO: complete it
			lattice = LatticeGen("HCP");
			lattice.GenerateC(c);
		}
		else {
			std::cout << "Undefined networks:\n";
		}


	}

	double NumDensity = polyC.NumParticle() / polyC.PeriodicVolume();
	if (Lengthscale.compare("D") == 0) {
		factor = D / polyC.NearestParticleDistance(0);
	}
	else if (Lengthscale.compare("rho") == 0) {
		factor = pow(rho / NumDensity, 1.0 / polyC.GetDimension());
	}
	else if (Lengthscale.compare("a") == 0) {
		factor = 1.0 / sqrt(polyC.GetBasisVector(0).Modulus2());
	}
	polyC.Rescale(factor);
	filename = NetworkName;

}


void GetPeriodicNetworks(std::istream & ifile, std::ostream & ofile, Dispersion &polyC, std::string & filename) {
	double	rho = 0,	/**< number density */
		D = 0,			/**< nearest-neighbor distances*/
		a = 0;			/**< lattice constant */

	std::string NetworkName, Lengthscale;

	ofile << "Network name = ";	ifile >> NetworkName;
	ofile << "Lengthscale (D/ rho/ a) = "; ifile >> Lengthscale;

	if (Lengthscale.compare("D") == 0) {
		ofile << "Nearest Neighbor distance, D = ";	ifile >> D;
	}
	else if (Lengthscale.compare("rho") == 0) {
		ofile << "Number density, rho = ";	ifile >> rho;
	}
	else if (Lengthscale.compare("a") == 0) {
		ofile << "Lattice constant = ";	ifile >> a;
	}
	else {
		ofile << "Choose a correct length scale !\n";
		return ;
	}

	LatticeGen lattice;
	Configuration c;
	double factor = 1.0;

	/*Define periodic network*/
	if (NetworkName.compare("custom") == 0) {
		ofile << "custom: is not implemented yet.\n";
		//TODO: input for custom periodic networks.
	}
	else {
		/*=========2D========*/
		if (NetworkName.compare("SquareNetwork") == 0) {
			lattice = LatticeGen("Square");
			lattice.GenerateC(c);
			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1);
			std::vector<GeometryVector> vertices;
			vertices.emplace_back(0.5*(v1 + v2));
			vertices.emplace_back(0.5*((-1.0)*v1 + v2));
			vertices.emplace_back(0.5*((-1.0)*v1 + (-1.0)*v2));
			vertices.emplace_back(0.5*(v1 - v2));

			polyC = Dispersion(c, Polygon(vertices).clone());
		}
		else if (NetworkName.compare("OctagonalNetwork") == 0) {
			lattice = LatticeGen("Square");
			lattice.GenerateC(c);
			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1);
			std::vector<GeometryVector> vertices;
			double f = 1.0 / (2.0 + sqrt(2.0));
			/*The first particle (square)*/
			vertices.emplace_back(f*v1);
			vertices.emplace_back(f*v2);
			vertices.emplace_back(-f * v1);
			vertices.emplace_back(-f * v2);
			polyC = Dispersion(c, Polygon(vertices).clone());
			polyC.MoveParticle(0, GeometryVector(0.0, 0.0));

			/*The second particle (octagon)*/
			GeometryVector x0 = 0.5*(v1 + v2);
			vertices[0] = f * v1 - x0;
			vertices[1] = (1.0 - f)*v1 - x0;
			vertices[2] = v1 + f * v2 - x0;
			vertices[3] = v1 + (1.0 - f)*v2 - x0;
			vertices.emplace_back((1.0 - f)*v1 + v2 - x0);
			vertices.emplace_back(f*v1 + v2 - x0);
			vertices.emplace_back((1.0 - f)*v2 - x0);
			vertices.emplace_back(f*v2 - x0);
			polyC.Insert(Polygon(vertices).clone(), GeometryVector(0.5, 0.5));

		}
		else if (NetworkName.compare("RhombicNetwork") == 0) {
			lattice = LatticeGen("Triangular");
			lattice.GenerateC(c);
			/*The first particle*/
			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1),
				x0 = 0.5*(v1+v2);
			std::vector<GeometryVector> vertices;
			vertices.emplace_back(0.0*v1 - x0);
			vertices.emplace_back(v1 - x0);
			vertices.emplace_back(v1 + v2 - x0);
			vertices.emplace_back(v2 - x0);
			polyC = Dispersion(c, Polygon(vertices).clone());
			polyC.MoveParticle(0, GeometryVector(0.5, 0.5));
		}
		else if (NetworkName.compare("TriangularNetwork") == 0) {
			lattice = LatticeGen("Triangular");
			lattice.GenerateC(c);
			/*The first particle*/
			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1);
			std::vector<GeometryVector> vertices;
			vertices.emplace_back(v1);
			vertices.emplace_back(v2);
			vertices.emplace_back(0.0*v1);
			polyC = Dispersion(c, Polygon(vertices).clone());
			polyC.MoveParticle(0, GeometryVector(1.0 / 3.0, 1.0 / 3.0));

			/*The second particle*/
			vertices[0] = v1;
			vertices[1] = v1 + v2;
			vertices[2] = v2;
			polyC.Insert(Polygon(vertices).clone(), GeometryVector(2.0 / 3.0, 2.0 / 3.0));
		}
		else if (NetworkName.compare("HoneycombNetwork") == 0) {
			lattice = LatticeGen("Triangular");
			lattice.GenerateC(c);

			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1);
			std::vector<GeometryVector> vertices;
			vertices.emplace_back(1.0 / 3.0*(v1 + v2));
			vertices.emplace_back(1.0 / 3.0*((-1)*v1 + 2 * v2));
			vertices.emplace_back(1.0 / 3.0*((-2)*v1 + 1 * v2));
			vertices.emplace_back(1.0 / 3.0*((-1)*v1 - 1 * v2));
			vertices.emplace_back(1.0 / 3.0*((1)*v1 - 2 * v2));
			vertices.emplace_back(1.0 / 3.0*((2)*v1 - 1 * v2));

			polyC = Dispersion(c, Polygon(vertices).clone());
		}
		else if (NetworkName.compare("KagomeNetwork") == 0) {
			lattice = LatticeGen("Triangular");
			lattice.GenerateC(c);
			/*The first particle (triangle)*/
			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1);
			std::vector<GeometryVector> vertices;
			GeometryVector x0 = 1.0 / 6.0*(v1 + v2); //center of mass...
			vertices.emplace_back(0.5*v1 - x0);
			vertices.emplace_back(0.5*v2 - x0);
			vertices.emplace_back(0.0*v1 - x0);
			polyC = Dispersion(c, Polygon(vertices).clone());
			polyC.MoveParticle(0, GeometryVector(1.0 / 6.0, 1.0 / 6.0));

			/*The second particle (triangle) */
			x0 = 5.0 / 6.0*(v1 + v2);
			vertices[0] = v1 + 0.5*v2 - x0;
			vertices[1] = v1 + v2 - x0;
			vertices[2] = 0.5*v1 + v2 -x0;
			polyC.Insert(Polygon(vertices).clone(), GeometryVector(5.0 / 6.0, 5.0 / 6.0));
			
			/*The third particle (hexagon) */
			x0 = 0.5*(v1 + v2);
			vertices[0] = 0.5*v1 - x0;
			vertices[1] = v1 - x0;
			vertices[2] = v1 + 0.5*v2 - x0;
			vertices.emplace_back(v2 + 0.5*v1 - x0);
			vertices.emplace_back(v2 - x0);
			vertices.emplace_back(0.5*v2 - x0);
			polyC.Insert(Polygon(vertices).clone(), GeometryVector(0.5, 0.5));
		}
		/* ==========3D======== */
		else if (NetworkName.compare("CubicNetwork") == 0) {
			lattice = LatticeGen("Cubic");
			lattice.GenerateC(c);

			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1), v3 = c.GetBasisVector(2),
				x0 = 0.5*(v1+v2+v3);
			std::vector<GeometryVector> vertices;	std::vector<std::vector<size_t>> index;
			/*0*/	vertices.emplace_back(0.0*v1 - x0);
			/*1*/	vertices.emplace_back(v2 - x0);
			/*2*/	vertices.emplace_back(v1 + v2 - x0);
			/*3*/	vertices.emplace_back(v1 - x0);
			/*4*/	vertices.emplace_back(v3 - x0);
			/*5*/	vertices.emplace_back(v1 + v3 - x0);
			/*6*/	vertices.emplace_back(v2 + v3 - x0);
			/*7*/	vertices.emplace_back(v1 + v2 + v3 - x0);
			{
				index.emplace_back(std::vector<size_t> { 0, 1, 2, 3 });
				index.emplace_back(std::vector<size_t> { 0, 3, 5, 4 });
				index.emplace_back(std::vector<size_t> { 0, 4, 6, 1 });
				index.emplace_back(std::vector<size_t> { 1, 6, 7, 2 });
				index.emplace_back(std::vector<size_t> { 7, 5, 3, 2 });
				index.emplace_back(std::vector<size_t> { 4, 5, 7, 6 });
			}

			polyC = Dispersion(c, Polyhedron(vertices, index).clone());
		}
		else if (NetworkName.compare("FCCNetwork") == 0) {
			/* Particle = unit cell */
			lattice = LatticeGen("FCC");
			lattice.GenerateC(c);

			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1), v3 = c.GetBasisVector(2),
				x0 = 0.5*(v1+v2+v3);
			std::vector<GeometryVector> vertices;	std::vector<std::vector<size_t>> index;

			/*0*/	vertices.emplace_back(0.0*v1 - x0);
			/*1*/	vertices.emplace_back(v2 - x0);
			/*2*/	vertices.emplace_back(v1 + v2 - x0);
			/*3*/	vertices.emplace_back(v1 - x0);
			/*4*/	vertices.emplace_back(v3 - x0);
			/*5*/	vertices.emplace_back(v1 + v3 - x0);
			/*6*/	vertices.emplace_back(v2 + v3 - x0);
			/*7*/	vertices.emplace_back(v1 + v2 + v3 - x0);
			{
				index.emplace_back(std::vector<size_t> { 0, 1, 2, 3 });
				index.emplace_back(std::vector<size_t> { 0, 3, 5, 4 });
				index.emplace_back(std::vector<size_t> { 0, 4, 6, 1 });
				index.emplace_back(std::vector<size_t> { 1, 6, 7, 2 });
				index.emplace_back(std::vector<size_t> { 7, 5, 3, 2 });
				index.emplace_back(std::vector<size_t> { 4, 5, 7, 6 });
			}

			polyC = Dispersion(c, Polyhedron(vertices, index).clone());
		}
		else if (NetworkName.compare("BCCNetwork") == 0) {
			/* particle = unit cell */
			lattice = LatticeGen("BCC");
			lattice.GenerateC(c);

			GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1), v3 = c.GetBasisVector(2),
				x0 = 0.5*(v1+v2+v3);
			std::vector<GeometryVector> vertices;	std::vector<std::vector<size_t>> index;

			/*0*/	vertices.emplace_back(0.0*v1 - x0);
			/*1*/	vertices.emplace_back(v2 - x0);
			/*2*/	vertices.emplace_back(v1 + v2 - x0);
			/*3*/	vertices.emplace_back(v1 - x0);
			/*4*/	vertices.emplace_back(v3 - x0);
			/*5*/	vertices.emplace_back(v1 + v3 - x0);
			/*6*/	vertices.emplace_back(v2 + v3 - x0);
			/*7*/	vertices.emplace_back(v1 + v2 + v3 - x0);
			{
				index.emplace_back(std::vector<size_t> { 0, 1, 2, 3 });
				index.emplace_back(std::vector<size_t> { 0, 3, 5, 4 });
				index.emplace_back(std::vector<size_t> { 0, 4, 6, 1 });
				index.emplace_back(std::vector<size_t> { 1, 6, 7, 2 });
				index.emplace_back(std::vector<size_t> { 7, 5, 3, 2 });
				index.emplace_back(std::vector<size_t> { 4, 5, 7, 6 });
			}

			polyC = Dispersion(c, Polyhedron(vertices, index).clone());
		}
		else if (NetworkName.compare("HCPNetwork") == 0) {
			/* particle = unit cell */
			//TODO: complete it
			lattice = LatticeGen("HCP");
			lattice.GenerateC(c);
		}
		else if (NetworkName.compare("FCCNewtork_WS") == 0) {
			/* Particle = Wigner-Seitz cell */
			//TODO: complete
			//lattice = LatticeGen("FCC");
			//lattice.GenerateP(sphC);

			//GeometryVector v1 = c.GetBasisVector(0), v2 = c.GetBasisVector(1), v3 = c.GetBasisVector(2);
			//std::vector<GeometryVector> vertices;	std::vector<std::vector<size_t>> index;

		}
		else if (NetworkName.compare("BCCNetwork_WS") == 0) {
			//TODO: complete it
			lattice = LatticeGen("BCC");
			lattice.GenerateC(c);
		}
		else if (NetworkName.compare("HCPNetwork_WS") == 0) {
			//TODO: complete it
			lattice = LatticeGen("HCP");
			lattice.GenerateC(c);
		}
		else {
			ofile << "Undefined networks:\n";
		}


	}

	double NumDensity = polyC.NumParticle() / polyC.PeriodicVolume();
	if (Lengthscale.compare("D") == 0) {
		factor = D / polyC.NearestParticleDistance(0);
	}
	else if (Lengthscale.compare("rho") == 0) {
		factor = pow(rho / NumDensity, 1.0 / polyC.GetDimension());
	}
	else if (Lengthscale.compare("a") == 0) {
		factor = 1.0 / sqrt(polyC.GetBasisVector(0).Modulus2());
	}
	polyC.Rescale(factor);
	filename = NetworkName;
}



void PeriodicNetworkParameters::ComputeBV_Fourier(const Dispersion & network, double phi, std::vector<GeometryVector> & result) {
	if (!Init_Parameter_Fourier) {
		std::cerr << "PeriodicNetworkParameters is not initialized\n";
		return;
	}
		
	Dispersion config(network);
	config.Rescale_prtOnly(pow(phi/network.PackingFraction(), 1.0 / network.GetDimension()));
	int d = config.GetDimension();
	double rho_Lambda = 1.0 / config.PeriodicVolume();

	/** Initialize the wavevectors at which the surface coefficient is computed. */
	if (this->Ks.size() == 0) {
		this->Ks.reserve((int)std::round(HyperSphere_Volume(d, this->kMax / (2.0*pi)) / rho_Lambda));
		Ks = GetKs(config, kMax, kMax);	/* A half of reciprocal lattices vectors */
	}

	//Initialize results
	{
		result.resize(0);
		GeometryVector k(kMin, 0, 0, 0);
		for (double K = kMin; K < kMax; K += dK) {
			result.emplace_back(k);
			k = k + GeometryVector(dK, 0, 0, 0);
		}

		result.emplace_back(k);
	}
	omp_set_num_threads(this->NumThreads);
#pragma omp parallel for schedule(guided) num_threads(this->NumThreads)
	for (size_t i = 0; i < Ks.size(); i++) {
		GeometryVector k = Ks[i];
		double k_modulus = sqrt((k.Modulus2()));
		double term = SpectralDensity_Dispersion(config, k);
		for (int idx = 0; idx < d + 1; idx++)
			term /= k_modulus;
		//int idx = (int)floor(log(k_modulus / Kmin) / logBinSize);
		int idx = (int)ceil((k_modulus - kMin) / dK);;
		idx = (idx <= 0) ? 0 : idx;
		if ((size_t)idx < result.size()) {
#pragma omp atomic
			result[idx].x[1] += term;

		}
	}

	//Cummulative sum!
	for (size_t j = 1; j < result.size(); j++)
		result[j].x[1] += result[j - 1].x[1];


	double factor = 2.0* pow(2.0, d)*rho_Lambda*boost::math::tgamma(1.0 + d / 2.0)*boost::math::tgamma(1.0 + d / 2.0) / ::pi;/* / pow(D, d + 1);*/
	for (size_t j = 0; j < result.size(); j++)
		result[j].x[1] *= factor;

}

void PeriodicNetworkParameters::ComputeBV_Fourier(const SpherePacking & network, double phi, std::vector<GeometryVector> & result) {
	if (!Init_Parameter_Fourier) {
		std::cerr << "PeriodicNetworkParameters is not initialized\n";
		return;
	}

	SpherePacking config(network);
	config.Rescale_prtOnly(pow(phi/network.PackingFraction(), 1.0 / network.GetDimension()));
	int d = config.GetDimension();
	double rho_Lambda = 1.0 / config.PeriodicVolume();

	/** Initialize the wavevectors at which the surface coefficient is computed. */
	if (this->Ks.size() == 0) {
		this->Ks.reserve((int)std::round(HyperSphere_Volume(d, this->kMax / (2.0*pi)) / rho_Lambda));
		Ks = GetKs(config, kMax, kMax);	/* A half of reciprocal lattices vectors */
	}

	//Initialize results
	{
		result.resize(0);
		GeometryVector k(kMin, 0, 0, 0);
		for (double K = kMin; K < kMax; K += dK) {
			result.emplace_back(k);
			k = k + GeometryVector(dK, 0, 0, 0);
		}

		result.emplace_back(k);
	}
	omp_set_num_threads(this->NumThreads);
#pragma omp parallel for schedule(guided) num_threads(this->NumThreads)
	for (size_t i = 0; i < Ks.size(); i++) {
		GeometryVector k = Ks[i];
		double k_modulus = sqrt((k.Modulus2()));
		double term = SpectralDensity(config, k);
		for (int idx = 0; idx < d + 1; idx++)
			term /= k_modulus;
		//int idx = (int)floor(log(k_modulus / Kmin) / logBinSize);
		int idx = (int)ceil((k_modulus - kMin) / dK);;
		idx = (idx <= 0) ? 0 : idx;
		if ((size_t)idx < result.size()) {
#pragma omp atomic
			result[idx].x[1] += term;

		}
	}

	//Cummulative sum!
	for (size_t j = 1; j < result.size(); j++)
		result[j].x[1] += result[j - 1].x[1];


	double factor = 2.0* pow(2.0, d)*rho_Lambda *boost::math::tgamma(1.0 + d / 2.0)*boost::math::tgamma(1.0 + d / 2.0) / ::pi;/* / pow(D, d + 1);*/
	for (size_t j = 0; j < result.size(); j++)
		result[j].x[1] *= factor;

}

void PeriodicNetworkParameters::ComputeBV_Direct(const SpherePacking & network, double phi, std::vector<GeometryVector> & result, std::vector<GeometryVector> & LV) {
	if (!Init_Parameter_Direct) {
		std::cerr << "PeriodicNetworkParameters is not initialized\n";
		return;
	}

	SpherePacking config(network);
	config.Rescale_prtOnly(pow(phi/ network.PackingFraction(), 1.0 / network.GetDimension()));
	int d = config.GetDimension();

	/*Initialize results*/
	{
		result.resize(0);
		for (double R = Rmin; R < Rmax; R += intv_R) {
			result.emplace_back(R, 0,0,0);
		}
		result.emplace_back(Rmax, 0, 0, 0);
	}

	std::vector<GeometryVector> LVFV;
	auto getConfig = [&network](size_t i)->SpherePacking {return network; };
	MC_vfv(getConfig, 1, LVFV, this->Rmax, this->dR, this->num_sampling_windows, this->NumThreads);

	int idx = 0;
	double sum = 0.0L, error = 0.0L;
	double first_val = 0.0L;	/* cache the original (scaled) LVFV[0].x[1] for the trapezoid left endpoint */
	for (size_t i = 0; i < LVFV.size(); i++) {
		double factor = d * pow(LVFV[i].x[0], d + 1.0);
		LVFV[i].x[1] *= factor;
		LVFV[i].x[2] *= factor;
		LVFV[i].x[3] = 0.0L;
		if (i == 0)
			first_val = LVFV[0].x[1];
		sum += LVFV[i].x[1];
		if (LVFV[i].x[1] != 0.0L)
			error += LVFV[i].x[2] * LVFV[i].x[2] / (LVFV[i].x[1] * LVFV[i].x[1]);

		LVFV[i].x[1] = (sum - 0.5*(first_val + LVFV[i].x[1])) / (i + 1.0);
		LVFV[i].x[2] = sqrt(error / (i + 1.0))*LVFV[i].x[1];
		if ((size_t)idx < result.size() && result[idx].x[0] <= LVFV[i].x[0] && LVFV[i].x[0] < result[idx].x[0] + dR) {
			result[idx].x[1] = LVFV[i].x[1];
			result[idx].x[2] = LVFV[i].x[2];
			idx++;
		}
	}
	LV = std::vector<GeometryVector>(LVFV);
}

void PeriodicNetworkParameters::ComputeVV_Fourier(const Dispersion & network, double phi, std::vector<GeometryVector>& result)
{
	if (!Init_Parameter_Fourier) {
		std::cerr << "PeriodicNetworkParameters is not initialized\n";
		return;
	}
	if (!Init_Parameter_Direct) {
		std::cerr << "PeriodicNetworkParameters is not initialized\n";
		return;
	}

	Dispersion config(network);
	config.Rescale_prtOnly(pow(phi / network.PackingFraction(), 1.0 / network.GetDimension()));
	int d = config.GetDimension();
	double rho_Lambda = 1.0 / config.PeriodicVolume();

	/** Initialize the wavevectors at which the surface coefficient is computed. */
	if (this->Ks.size() == 0) {
		this->Ks.reserve((int)std::round(HyperSphere_Volume(d, this->kMax / (2.0*pi)) / rho_Lambda));
		Ks = GetKs(config, kMax, kMax);	/* A half of reciprocal lattices vectors */
	}

	//Initialize results
	{
		result.resize(0);
		GeometryVector temp(0, 0, 0, 0);
		for (double R = Rmin; R < Rmax; R += dR) {
			temp.x[0] = R;
			result.emplace_back(temp);
		}
	}

	std::vector<GeometryVector> chik;
	std::vector<GeometryVector> Xs(this->Ks.size(), GeometryVector(0.0, 0.0));
	omp_set_num_threads(this->NumThreads);
#pragma omp parallel for schedule(guided) num_threads(this->NumThreads)
	for (size_t i = 0; i < Ks.size(); i++) {
		GeometryVector k = Ks[i];
		double term = SpectralDensity_Dispersion(config, k);

		Xs[i].x[0] = sqrt(k.Modulus2());
		Xs[i].x[1] = term;
	}

	auto CompareFunc = [](const GeometryVector & left, const GeometryVector & right) ->bool
	{
		return left.x[0] < right.x[0];
	};
	std::sort(Xs.begin(), Xs.end(), CompareFunc);

	if (Ks.empty() || Xs.empty())
		return;
	chik.reserve(int(Ks.size() / HyperSphere_SurfaceArea(d, Ks[Ks.size() - 1].x[0])));
	int j = 0;
	chik.emplace_back(Xs[0]);
	double tracer = Xs[0].x[0];
	for (size_t i = 1; i < Xs.size(); i++) {
		/*wavenumber increased */
		if (Xs[i].x[0] - tracer > 1e-7) {
			j++;	tracer = Xs[i].x[0];
			chik.emplace_back(Xs[i].x[0], Xs[i].x[1]);
		}
		/*wavenumber is the same */
		else {
			chik[j].x[1] += Xs[i].x[1];
		}
	}
	//WriteFunction(chik, "test-square");
	Xs.clear();
	/* chik[j] = {wavenumber, \sum_{\abs{\vect{k}} = wavenumber} Ks[i].x[1]} */


	/*  --- compute the local volume-fraction variance. --- */
	double prefactor = 2.0* boost::math::tgamma(1.0 + d / 2.0)*boost::math::tgamma(1.0 + d / 2.0) / (config.PeriodicVolume()*pow(pi, d));
	for (size_t i = 0; i < result.size(); i++) {
		double R = result[i].x[0], 
			factor = prefactor / pow(R*R, d),
			sum = 0.0L;
		
#pragma omp parallel for 
		for (size_t j = 0; j < chik.size(); j++) {
			double m = Mtilde(chik[j].x[0], R, d),
				m2 = m * m;
#pragma omp atomic
			sum += m2 * chik[j].x[1];
		}
		
		/*sum *= factor;*/
		result[i].x[1] = sum * factor;
	}
}

void PeriodicNetworkParameters::ComputeVV_Fourier(const SpherePacking & network, double phi_void, std::vector<GeometryVector>& result)
{
	if (!Init_Parameter_Fourier) {
		std::cerr << "PeriodicNetworkParameters is not initialized\n";
		return;
	}
	if (!Init_Parameter_Direct) {
		std::cerr << "PeriodicNetworkParameters is not initialized\n";
		return;
	}

	SpherePacking config(network);
	config.Rescale_prtOnly(pow(phi_void / network.PackingFraction(), 1.0 / network.GetDimension()));
	int d = config.GetDimension();
	double rho_Lambda = 1.0 / config.PeriodicVolume();

	/** Initialize the wavevectors at which the surface coefficient is computed. */
	if (this->Ks.size() == 0) {
		this->Ks.reserve((int)std::round(HyperSphere_Volume(d, this->kMax / (2.0*pi)) / rho_Lambda));
		Ks = GetKs(config, kMax, kMax);	/* A half of reciprocal lattices vectors */
	}

	//Initialize results
	{
		result.resize(0);
		GeometryVector temp(0, 0, 0, 0);
		for (double R = Rmin; R < Rmax; R += dR) {
			temp.x[0] = R;
			result.emplace_back(temp);
		}
	}

	/* ---compute the spectral density in advance. ---*/
	std::vector<GeometryVector> chik;
	std::vector<GeometryVector> Xs(this->Ks.size(), GeometryVector(0.0,0.0));
	omp_set_num_threads(this->NumThreads);
#pragma omp parallel for schedule(guided) num_threads(this->NumThreads)
	for (size_t i = 0; i < Ks.size(); i++) {
		GeometryVector k = Ks[i];
		double term = SpectralDensity(config, k);
		
		Xs[i].x[0] = sqrt(k.Modulus2());
		Xs[i].x[1] = term;
	}

	auto CompareFunc = [](const GeometryVector & left, const GeometryVector & right) ->bool
	{
		return left.x[0] < right.x[0];
	};
	std::sort(Xs.begin(), Xs.end(), CompareFunc);

	if (Ks.empty() || Xs.empty())
		return;
	chik.reserve(int(Ks.size() / HyperSphere_SurfaceArea(d, Ks[Ks.size() - 1].x[0])));
	int j = 0;
	chik.emplace_back(Xs[0]);
	double tracer = Xs[0].x[0];
	for (size_t i = 1; i < Xs.size(); i++) {
		/*wavenumber increased */
		if (Xs[i].x[0] - tracer > 1e-7) {
			j++;	tracer = Xs[i].x[0];
			chik.emplace_back(Xs[i].x[0], Xs[i].x[1]);
		}
		/*wavenumber is the same */
		else {
			chik[j].x[1] += Xs[i].x[1];
		}
	}
	//WriteFunction(chik, "test-square");
	Xs.clear();
	/* chik[j] = {wavenumber, \sum_{\abs{\vect{k}} = wavenumber} Ks[i].x[1]} */


	/*  --- compute the local volume-fraction variance. --- */
	double prefactor = 2.0* boost::math::tgamma(1.0 + d / 2.0)*boost::math::tgamma(1.0 + d / 2.0) / (config.PeriodicVolume()*pow(pi, d));
	for (size_t i = 0; i < result.size(); i++) {
		double R = result[i].x[0],
			factor = prefactor/pow(R*R,d),
			sum = 0.0L;

#pragma omp parallel for 
		for (size_t j = 0; j < chik.size(); j++) {
			double m = Mtilde(chik[j].x[0], R, d),
				m2 = m * m;
#pragma omp atomic
			sum += m2 * chik[j].x[1];
		}

		/*sum *= factor;*/
		result[i].x[1] = sum * factor;
	}
}

void PeriodicNetworkParameters::ComputeVV_Direct(const SpherePacking & network, double phi_void, std::vector<GeometryVector>& result)
{
	if (!Init_Parameter_Direct) {
		std::cerr << "PeriodicNetworkParameters is not initialized\n";
		return;
	}

	//Initialize results
	{
		result.resize(0);
		GeometryVector temp(0, 0, 0, 0);
		for (double R = dR; R < Rmax; R += dR) {
			temp.x[0] = R;
			result.emplace_back(temp);
		}
	}

	int num_ensemble = 50,
		num_samp = this->num_sampling_windows/ num_ensemble;

	auto getc = [&network, &phi_void](size_t i)->SpherePacking{
		double factor = phi_void / network.PackingFraction();
		factor = pow(factor, 1.0 / (double)network.GetDimension());
		SpherePacking c(network);
		c.Rescale_prtOnly(factor);
		return c;
	};

#pragma omp parallel for schedule(guided) num_threads(this->NumThreads)
	for (int i = 0; i < num_ensemble; i++) {
		std::vector<GeometryVector> result_t;
		MC_vfv(getc, 1, result_t, this->Rmax, this->dR, num_samp, 1);
		
		for (size_t j = 0; j < result_t.size(); j++) {
#pragma omp atomic
			result[j].x[1] += result_t[j].x[1];
#pragma omp atomic
			result[j].x[2] += result_t[j].x[1] * result_t[j].x[1];
		}
	}

	for (size_t i = 0; i < result.size(); i++) {
		result[i].x[1] /= (double)num_ensemble;
		result[i].x[2] /= (double)num_ensemble;
		result[i].x[2] = sqrt(abs(result[i].x[2] - result[i].x[1]* result[i].x[1])/(num_ensemble-1.0));
	}
}


/**	 @see PerioidcNetwork::CheckFundamentalCell()
 *	 It doesn't check whether particles are overlapping or not.
 */
bool PeriodicNetwork::CheckFundamentalCell() {
	bool result = true;

	// the number of particles
	result &= (this->c.NumParticle()) == (this->Particles.size());
	// volume fraction
	{
		double VolumeFraction = 0.0;
		for (size_t i = 0; i < this->Particles.size(); i++)
			VolumeFraction += this->Particles[i]->GetVolume();

		VolumeFraction /= this->c.PeriodicVolume();
		if (this->phi < 0) {//Uninitialized volume fraction => initialize volume fraction.
			this->phi = VolumeFraction;
		}

		result &= (this->phi == VolumeFraction);
	}
	return result;
}


/**	Resize particles.
 *	@param volumefraction	Resize particles in a target packing fraction 
 */
void PeriodicNetwork::SetPhi(double volumefraction) {
	double factor = pow(volumefraction / this->phi, 1.0 / this->c.GetDimension());
	for (size_t iter = 0; iter < this->Particles.size(); iter++) {
		this->Particles[iter]->Scale(factor);
	}
	this->phi = volumefraction;
}

/** Print particles.
 *	@param ofile	An output stream.
 */
void PeriodicNetwork::Print(std::ostream & ofile) const {
	int d = this->c.GetDimension();
	ofile << "Basis: " << d << "\n";
	for(int i=0; i<d; i++)
		ofile << this->c.GetBasisVector(i);

	ofile << "\n";
	for (size_t i = 0; i < this->Particles.size(); i++) {
		ofile << "x_" << i << ":\t"<<this->c.GetCartesianCoordinates(i);
		ofile << "Particle "<<i<<":\n";
		this->Particles[i]->Print(ofile);
	}
}

/** 
 */
void  PeriodicNetwork::ComputeBV_Fourier(double Kmin, double Kmax, double dk, std::vector<GeometryVector> & result, std::vector<GeometryVector> & Ks) {
	double	rho_Lambda = 1.0 / this->c.PeriodicVolume(),	/**< number density of the "Bravais lattice points" */
			D = this->c.NearestParticleDistance(0);						/**< nearest-neighbor distances*/

	int d = this->c.GetDimension();
	
	//std::vector<GeometryVector> Ks;	
	if (Ks.size() == 0) {
		Ks.reserve((int)std::round(HyperSphere_Volume(d, Kmax / (2.0*pi)) / rho_Lambda));
		Ks = GetKs(this->c, Kmax, Kmax);	/* A half of reciprocal lattices vectors */
	}
	
	//Initialize results
	{	
		result.resize(0);
		GeometryVector k(Kmin, 0, 0, 0);
		for (double K = Kmin; K < Kmax; K += dk) {
			result.emplace_back(k);
			k = k + GeometryVector(dk, 0, 0, 0);
		}
		
		result.emplace_back(k);
	}
	omp_set_num_threads(this->num_threads);
#pragma omp parallel for schedule(guided) num_threads(this->num_threads)
	for (size_t i = 0; i < Ks.size(); i++) {
		GeometryVector k = Ks[i];
		double k_modulus = sqrt((k.Modulus2()));
		double term = GetSquaredFormFactor(k);
		for (int idx = 0; idx < d + 1; idx++)
			term /= k_modulus;
		//int idx = (int)floor(log(k_modulus / Kmin) / logBinSize);
					int idx = (int)ceil((k_modulus - Kmin) / dk);;
		idx = (idx <= 0) ? 0 : idx;
					if ((size_t)idx < result.size()) {
#pragma omp atomic
			result[idx].x[1] += term;

		}
	}

		//Cummulative sum!
		for (size_t j = 1; j < result.size(); j++)
			result[j].x[1] += result[j - 1].x[1];


		double factor = pow(2.0, d)*rho_Lambda*rho_Lambda*boost::math::tgamma(1.0 + d / 2.0)*boost::math::tgamma(1.0 + d / 2.0) / ::pi / pow(D, d + 1);
		factor *= 2;	//Because ks do not include parity pairs.
		for (size_t j = 0; j < result.size(); j++)
			result[j].x[1] *= factor;

}



/**	A CLI interface to compute BV for some periodic networks.
 *
 * 2D: SquareNetwork,	TriangularNetwork,	HexagonalNetwork, 
 * 3D: CubicNetwork,	BCC,	FCC,
 */
void ComputeBV_CLI() {
	std::ostream & ofile = std::cout;
	std::istream & ifile = std::cin;

	double	phi_init,	/**< initial volume fraction */
			phi_fin,	/**< final volume fraction*/
			d_phi;		/**< change in volume fraction */
	bool DoDirectComputation = false;

	/* Parameters for computing B_V */
	PeriodicNetworkParameters Param;
	
	std::string ParticleType;
	ofile << "Particle Type (sphere/polygon) = "; ifile >> ParticleType;
	ofile << "Do compute Bv in the direct space? (sphere only) = "; ifile >> DoDirectComputation;
	

	if (DoDirectComputation) {
		if (ParticleType.compare("sphere") != 0) {
			DoDirectComputation = false;
			ofile << "Direct computation for polygonal particles is not implemented\n"
				<< "We compute Bv in the Fourier space\n";
		}
		else {
			double Rmax = 400, dR = 0.1, Rmin = 20, intv_R = 20;
			ofile << "Maximal radius = ";	ifile >> Rmax;
			ofile << "dR = ";				ifile >> dR;
			ofile << "Minimal radius for average = ";	ifile >> Rmin;
			ofile << "Interval of R for average = ";	ifile >> intv_R;

			Param.SetParam_Direct(Rmax, dR, Rmin, intv_R);
		}
	}
	if (!DoDirectComputation){
		double kmin = 500, kmax = 1000, dk = 100;
		ofile << "Maximal wavenumber = ";	ifile >> kmax;
		ofile << "Minimal wavenumber for summation = ";	ifile >> kmin;
		ofile << "Change in wavenumber = ";	ifile >> dk;

		Param.SetParam_Fourier(kmax, kmin, dk);
	}
	ofile << "The number of threads = "; ifile >> Param.NumThreads;

	ofile << "Initial volume fraction = ";	ifile >> phi_init;
	ofile << "Final volume fraction = ";	ifile >> phi_fin;
	ofile << "Change in volume fraction = ";	ifile >> d_phi;

	char name[300];
	

	if (ParticleType.compare("sphere") == 0) {
		/* Determine the fundamental cell */
		std::string keepgoing = "y";

		for (; keepgoing.compare("y") == 0 || keepgoing.compare("Y") == 0;) {
			std::fstream save;
			std::string fileprefix;

			SpherePacking sphC; GetPeriodicPackings(ifile, ofile, sphC, fileprefix);
			phi_fin = sphC.PackingFraction();
			ofile << "Computing: sphere packing-" << fileprefix << "\n";
			if (DoDirectComputation) {
				ofile << "\tin direct space, ";
				/*Compute B_V from local volume-fraction variance */
				std::sprintf(name, "%s_Bv_Direct.txt", fileprefix.c_str());
				save.open(name, std::fstream::out);
				save << "0";
				for (double R = Param.Rmin; R <= Param.Rmax; R += Param.intv_R)
					save << "\t" << R;
				save << "\n";

				for (double phi = phi_init; phi < phi_fin + d_phi; phi += d_phi) {
					phi = std::min(phi_fin, phi);
					SpherePacking config(sphC);	config.Rescale_prtOnly(pow(phi / phi_fin, 1.0 / sphC.GetDimension()));
					std::vector<GeometryVector> result, LV;
					Param.ComputeBV_Direct(config, phi, result, LV);
					WriteFunction(LV, (fileprefix + "_runningav_"+std::to_string(phi)).c_str());
					save << phi;
					for (size_t i = 0; i < result.size(); i++)
						save << "\t" << result[i].x[1];
					save << "\n";
				}
				ofile << "done\n";
			}
			else {
				ofile << "\tin Fourier space, ";
				/*Compute B_V from the Fourier sum */
				std::sprintf(name, "%s_Bv_Fourier.txt", fileprefix.c_str());
				save.open(name, std::fstream::out);
				save << "0";
				for (double k = Param.kMin; k <= Param.kMax; k += Param.dK)
					save << "\t" << k;
				save << "\n";

				for (double phi = phi_init; phi < phi_fin + d_phi; phi += d_phi) {
					phi = std::min(phi_fin, phi);
					SpherePacking config(sphC);	config.Rescale_prtOnly(pow(phi / phi_fin, 1.0 / sphC.GetDimension()));
					std::vector<GeometryVector> result;
					Param.ComputeBV_Fourier(config, phi, result);
					save << phi;
					for (size_t i = 0; i < result.size(); i++)
						save << "\t" << result[i].x[1];
					save << "\n";
				}
				ofile << "done\n";
				Param.Ks.clear();/*Iinitialize the list of wavevectors.*/
			}
			ofile << "Compute other sphere packings? (y/n)\n";
			ifile >> keepgoing;
			save.close();
			ofile << "\n\n\n";
		}
	}
	else if( ParticleType.compare("polygon") == 0 ) {
		/* Determine the fundamental cell */
		std::string fileprefix;
		std::string keepgoing = "y";

		for (; keepgoing.compare("y") == 0 || keepgoing.compare("Y") == 0;) {
			std::fstream save;
			std::string fileprefix;

			Dispersion C; GetPeriodicNetworks(ifile, ofile, C, fileprefix);
			ofile << "Computing: periodic network-"<<fileprefix<<"\n";

			if (DoDirectComputation) {
				//TODO: !!
			}
			else {
				/*Compute B_V from the Fourier sum */
				ofile << "\tin Fourier space, ";
				std::sprintf(name, "%s_Bv_Fourier.txt", fileprefix.c_str());
				save.open(name, std::fstream::out);
				save << "0";
				for (double k = Param.kMin; k <= Param.kMax; k += Param.dK)
					save << "\t" << k;
				save << "\n";

				for (double phi = phi_init; phi < phi_fin + d_phi; phi += d_phi) {
					phi = std::min(phi_fin, phi);
					Dispersion config(C);	config.Rescale_prtOnly(pow(phi / C.PackingFraction(), 1.0 / C.GetDimension()));
					std::vector<GeometryVector> result;
					Param.ComputeBV_Fourier(config, phi, result);
					save << phi;
					for (size_t i = 0; i < result.size(); i++)
						save << "\t" << result[i].x[1];
					save << "\n";
				}
				ofile << "done\n";
				Param.Ks.clear();/*Iinitialize the list of wavevectors.*/
			}
			ofile << "Compute other periodic networks? (y/n)\n";
			ifile >> keepgoing;
			ofile << "\n\n\n";
			save.close();
		}
	}
	else {
		std::cerr << "Wrong particle type\n";
		return;
	}
	

}