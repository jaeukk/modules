#ifndef PHONONFREQUENCY_INCLUDED
#define PHONONFREQUENCY_INCLUDED

#include <GeometryVector.h>
#include <Potential.h>
#include <PeriodicCellList.h>
#include <vector>
/** \brief Compute the squared phonon frequencies of a crystal at a wave vector. (CLAUDE)
 *  Builds the dynamical matrix from the pair potential's second derivatives summed
 *  over neighbors, diagonalizes it, and returns the sorted squared frequencies.
 *  \param k        Wave vector at which the dynamical matrix is evaluated.
 *  \param crystal  Periodic crystal configuration (atom positions and lattice).
 *  \param pot      Pair potential providing the force constants (second derivatives).
 *  \param result2  Output, filled with the squared frequencies omega^2 in ascending order.
 *  \param Mass     Atomic mass used to scale the eigenvalues into omega^2. */
void PhononFrequency(const GeometryVector & k, const Configuration & crystal, PairPotential & pot, std::vector<double> & result2, double Mass=1.0);

/** \brief Compute and plot the phonon dispersion along a path of k-space end points. (CLAUDE)
 *  \param kEndPoints       Successive k-space points defining the band-structure path.
 *  \param crystal          Periodic crystal configuration.
 *  \param pot              Pair potential providing the force constants.
 *  \param Prefix           Filename prefix for the output data and figure files.
 *  \param Labels           High-symmetry point labels, one per k end point.
 *  \param Mass             Atomic mass passed to the frequency calculation.
 *  \param NumKPointsPerLine Number of sampled k points per segment of the path. */
void PlotPhononFrequencies(const std::vector<GeometryVector> & kEndPoints, const Configuration & crystal, PairPotential & pot, char * Prefix, std::vector<const char *> & Labels, double Mass=1.0, size_t NumKPointsPerLine=30);


extern double LowestFrequencySquared;//use it to broadcast calculation results to other functions

/** \brief Plot the phonon dispersion of a 3D hexagonal crystal along its standard path. (CLAUDE)
 *  \param crystal          Periodic crystal configuration (must be 3D).
 *  \param pot              Pair potential providing the force constants.
 *  \param Prefix           Filename prefix for the output data and figure files.
 *  \param Mass             Atomic mass passed to the frequency calculation.
 *  \param NumKPointsPerLine Number of sampled k points per segment of the path. */
void PlotPhononFrequencies_3DHexagonal(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=30);

/** \brief Plot the phonon dispersion of a 3D tetragonal crystal along its standard path. (CLAUDE)
 *  \param crystal          Periodic crystal configuration (must be 3D).
 *  \param pot              Pair potential providing the force constants.
 *  \param Prefix           Filename prefix for the output data and figure files.
 *  \param Mass             Atomic mass passed to the frequency calculation.
 *  \param NumKPointsPerLine Number of sampled k points per segment of the path. */
void PlotPhononFrequencies_3DTetragonal(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=30);

/** \brief Plot the phonon dispersion of a 3D cubic crystal along its standard path. (CLAUDE)
 *  \param crystal          Periodic crystal configuration (must be 3D).
 *  \param pot              Pair potential providing the force constants.
 *  \param Prefix           Filename prefix for the output data and figure files.
 *  \param Mass             Atomic mass passed to the frequency calculation.
 *  \param NumKPointsPerLine Number of sampled k points per segment of the path. */
void PlotPhononFrequencies_3DCubic(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=30);

/** \brief Plot the phonon dispersion of a 2D hexagonal crystal along its standard path. (CLAUDE)
 *  \param crystal          Periodic crystal configuration (must be 2D).
 *  \param pot              Pair potential providing the force constants.
 *  \param Prefix           Filename prefix for the output data and figure files.
 *  \param Mass             Atomic mass passed to the frequency calculation.
 *  \param NumKPointsPerLine Number of sampled k points per segment of the path. */
void PlotPhononFrequencies_2DHexagonal(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=90);

/** \brief Plot the phonon dispersion of a 2D square (cubic) crystal along its standard path. (CLAUDE)
 *  \param crystal          Periodic crystal configuration (must be 2D).
 *  \param pot              Pair potential providing the force constants.
 *  \param Prefix           Filename prefix for the output data and figure files.
 *  \param Mass             Atomic mass passed to the frequency calculation.
 *  \param NumKPointsPerLine Number of sampled k points per segment of the path. */
void PlotPhononFrequencies_2DCubic(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=90);

/** \brief Plot the phonon dispersion of a 2D rectangular crystal along its standard path. (CLAUDE)
 *  \param crystal          Periodic crystal configuration (must be 2D).
 *  \param pot              Pair potential providing the force constants.
 *  \param Prefix           Filename prefix for the output data and figure files.
 *  \param Mass             Atomic mass passed to the frequency calculation.
 *  \param NumKPointsPerLine Number of sampled k points per segment of the path. */
void PlotPhononFrequencies_2DRectangle(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=60);

/** \brief Plot the phonon dispersion of a 3D FCC crystal along its standard path. (CLAUDE)
 *  \param crystal          Periodic crystal configuration (must be 3D).
 *  \param pot              Pair potential providing the force constants.
 *  \param Prefix           Filename prefix for the output data and figure files.
 *  \param Mass             Atomic mass passed to the frequency calculation.
 *  \param NumKPointsPerLine Number of sampled k points per segment of the path. */
void PlotPhononFrequencies_3DFCC(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=20);

#endif