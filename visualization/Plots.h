#ifndef PLOTS_INCLUDED
#define PLOTS_INCLUDED

#include "ParameterSet.h"

#if defined(LATTICESUMSOLVER_INCLUDED)
#include "LatticeSumSolver.h"
/** \brief Plot the pair potential and lattice sums of candidate structures for a potential. (CLAUDE)
 *  \param pPot      Isotropic pair potential to plot.
 *  \param Solvers   Candidate-structure lattice-sum solvers to evaluate and plot.
 *  \param Prefix    Output file name prefix.
 *  \param Enphasis  Index of the solver curve to emphasize. */
void Plots(IsotropicPairPotential * pPot, const std::vector<LatticeSumSolver *> * Solvers, const char * Prefix, size_t Enphasis = 0);
/** \brief Plot the pair potential and lattice sums for a parameter set, writing its parameters too. (CLAUDE)
 *  \param Param     Parameter set describing the potential.
 *  \param Solvers   Candidate-structure lattice-sum solvers to evaluate and plot.
 *  \param Prefix    Output file name prefix.
 *  \param Enphasis  Index of the solver curve to emphasize. */
void Plots(const ParameterSet & Param, const std::vector<LatticeSumSolver *> * Solvers, const char * Prefix, size_t Enphasis = 0);
#endif

//parameter ParticleSizeRescale only works for 2D configurations
/** \brief Plot a 1D, 2D, or 3D configuration with its basis cell and particles. (CLAUDE)
 *  \param prefix               Output file name prefix.
 *  \param c                    Configuration to plot.
 *  \param remark               Optional caption drawn on the figure.
 *  \param ParticleSizeRescale  Rescaling factor for the drawn particle size (2D only). */
void Plot(std::string prefix, const Configuration & c, const std::string & remark="", double ParticleSizeRescale=1.0);


//require that all configurations in this vector have the same basis vectors
/** \brief Plot the directional (2D) structure factor averaged over several configurations. (CLAUDE)
 *  \param vc            Configurations to average; all must share the same basis vectors.
 *  \param limit         Number of reciprocal lattice steps to scan in each direction.
 *  \param prefix        Output file name prefix.
 *  \param LogScale      If true, plot on a log10 scale.
 *  \param MinSOverride  Lower clamp for S(k); 0 means no clamp.
 *  \param MaxSOverride  Upper clamp for S(k); 0 means no clamp.
 *  \param TileSize      Tile size used for the rendered colour map. */
void PlotDirectionalStructureFactor2D(const std::vector<Configuration> & vc, long limit, std::string prefix, bool LogScale, double MinSOverride, double MaxSOverride, double TileSize);
/** \brief Plot the directional (2D) structure factor of a single configuration. (CLAUDE)
 *  \param c             Configuration to analyze.
 *  \param limit         Number of reciprocal lattice steps to scan in each direction.
 *  \param prefix        Output file name prefix.
 *  \param LogScale      If true, plot on a log10 scale.
 *  \param MinSOverride  Lower clamp for S(k); 0 means no clamp.
 *  \param MaxSOverride  Upper clamp for S(k); 0 means no clamp.
 *  \param TileSize      Tile size used for the rendered colour map. */
void PlotDirectionalStructureFactor2D(const Configuration & c, long limit, std::string prefix, bool LogScale=true, double MinSOverride=0.0, double MaxSOverride=0.0, double TileSize=1.0);
/** \brief Plot the directional (2D) structure factor of a configuration (prefix-first overload). (CLAUDE)
 *  \param prefix        Output file name prefix.
 *  \param c             Configuration to analyze.
 *  \param limit         Number of reciprocal lattice steps to scan in each direction.
 *  \param LogScale      If true, plot on a log10 scale.
 *  \param MinSOverride  Lower clamp for S(k); 0 means no clamp.
 *  \param MaxSOverride  Upper clamp for S(k); 0 means no clamp.
 *  \param TileSize      Tile size used for the rendered colour map. */
inline void PlotDirectionalStructureFactor2D(std::string prefix, const Configuration & c, long limit, bool LogScale=true, double MinSOverride=0.0, double MaxSOverride=0.0, double TileSize=1.0)
{
	PlotDirectionalStructureFactor2D(c, limit, prefix, LogScale, MinSOverride, MaxSOverride, TileSize);
}

//Plot a path in a 2D configuration
/** \brief Plot a 2D configuration together with a path connecting selected particles. (CLAUDE)
 *  \param prefix  Output file name prefix.
 *  \param List    2D configuration to plot.
 *  \param remark  Optional caption drawn on the figure.
 *  \param path    Ordered list of particle indices forming the path. */
void PlotPath(std::string prefix, const Configuration & List, const std::string & remark, const std::vector<size_t> & path);

/** \brief Plot a 2D configuration with a displacement (arrow) field on each particle. (CLAUDE)
 *  \param prefix        Output file name prefix.
 *  \param List          2D configuration to plot.
 *  \param remark        Optional caption drawn on the figure.
 *  \param displacement  Flat (x, y) displacement components, two per particle. */
void PlotDisplacement(std::string prefix, const Configuration & List, const std::string & remark, const std::vector<double> & displacement);
/** \brief Plot a 2D periodic cell list as a sphere packing of uniform radius. (CLAUDE)
 *  \tparam add     Per-particle characteristics type of the cell list.
 *  \param prefix   Output file name prefix.
 *  \param List     Periodic cell list to plot.
 *  \param remark   Optional caption drawn on the figure.
 *  \param r        Radius assigned to every particle. */
template<typename add>
void PlotPacking2D(std::string prefix, const PeriodicCellList<add> & List, const std::string & remark, double r);
/** \brief Plot a 2D sphere packing with per-particle radii. (CLAUDE)
 *  \param prefix       Output file name prefix.
 *  \param List         Sphere packing to plot.
 *  \param sphereColor  MathGL colour code for the spheres.
 *  \param remark       Optional caption drawn on the figure. */
void PlotPacking2D(std::string prefix, const SpherePacking & List, const char * sphereColor = "b", const std::string & remark="");


/** \brief Plot a function (with error bars) to a bitmap/EPS file using MathGL. (CLAUDE)
 *  \param result         Data points; each x packs (x, y, ex, ey).
 *  \param OutFilePrefix  Output file name prefix.
 *  \param xLabel         X axis label.
 *  \param yLabel         Y axis label. */
void PlotFunction_MathGL(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::string & xLabel, const std::string & yLabel);

/** \brief Write one or more data series to a Grace (.agr) project file, with full axis control. (CLAUDE)
 *  \param presult        Array of NumSet data series.
 *  \param NumSet         Number of series in presult.
 *  \param OutFilePrefix  Output file name prefix (".agr" is appended).
 *  \param xLabel         X axis label.
 *  \param yLabel         Y axis label.
 *  \param legends        Legend text per series.
 *  \param Title          Plot title.
 *  \param MinX           Minimum X of the plot world.
 *  \param MaxX           Maximum X of the plot world.
 *  \param TickX          Major tick spacing on X.
 *  \param MinY           Minimum Y of the plot world.
 *  \param MaxY           Maximum Y of the plot world.
 *  \param TickY          Major tick spacing on Y. */
void PlotFunction_Grace(const std::vector<GeometryVector> * presult, size_t NumSet, const std::string & OutFilePrefix, const std::string & xLabel, const std::string & yLabel, const std::vector<std::string> & legends, const std::string & Title, double MinX, double MaxX, double TickX, double MinY, double MaxY, double TickY);
/** \brief Write one or more data series to a Grace (.agr) project file with auto-computed axis ranges. (CLAUDE)
 *  \param presult        Array of NumSet data series.
 *  \param NumSet         Number of series in presult.
 *  \param OutFilePrefix  Output file name prefix (".agr" is appended).
 *  \param xLabel         X axis label.
 *  \param yLabel         Y axis label.
 *  \param legends        Legend text per series.
 *  \param Title          Plot title. */
void PlotFunction_Grace(const std::vector<GeometryVector> * presult, size_t NumSet, const std::string & OutFilePrefix, const std::string & xLabel, const std::string & yLabel, const std::vector<std::string> & legends, const std::string & Title);
/** \brief Write a single data series to a Grace (.agr) project file. (CLAUDE)
 *  \param result         Data series to write.
 *  \param OutFilePrefix  Output file name prefix (".agr" is appended).
 *  \param xLabel         X axis label.
 *  \param yLabel         Y axis label.
 *  \param Title          Plot title. */
void PlotFunction_Grace(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix="temp", const std::string & xLabel="", const std::string & yLabel="", const std::string & Title="");

/** \brief Read data series from a Grace (.agr) file given its name prefix. (CLAUDE)
 *  \param result        Destination vector of series, filled from the file.
 *  \param InFilePrefix  Input file name prefix (".agr" is appended). */
void ReadGraceData(std::vector<std::vector<GeometryVector> > & result, const std::string & InFilePrefix);
/** \brief Read data series from an open Grace (.agr) input stream. (CLAUDE)
 *  \param result  Destination vector of series, filled from the stream.
 *  \param ifile   Input stream positioned at the start of the Grace file. */
void ReadGraceData(std::vector<std::vector<GeometryVector> > & result, std::istream & ifile);
#endif