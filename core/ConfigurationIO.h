#ifndef ConfigurationIO_INCLUDED
#define ConfigurationIO_INCLUDED

/** \file ConfigurationIO.h
 *  \brief File-format readers/writers for Configuration and SpherePacking
 *         (.pos, binary .Config, HOOMD XML, stealth, GROMACS .gro, raw
 *         coordinates, Torquato-Jiao packings, and plain .txt).  Split out of
 *         PeriodicCellList.h for readability. */

#include "PeriodicCellList.h"   // Configuration, PeriodicCellList<>
#include "SpherePacking.h"      // SpherePacking (TJ I/O + WriteConfiguration overloads)

/** \brief Write a Configuration to a .pos text file with the given prefix. 
 *  \param prefix File name prefix (".pos" is appended).
 *  \param List   Configuration to write. */
void Output(std::string prefix, const Configuration & List);

/** \brief Write a Configuration in .pos text format to an output stream. 
 *  \param out  Output stream to write to.
 *  \param List Configuration to write. */
void Output(std::ostream & out, const Configuration & List);

/** \brief Write a Configuration in binary format to an output stream. 
 *  \param out  Output stream to write to.
 *  \param List Configuration to write. */
inline void WriteBinary(std::ostream & out, const Configuration & List)
{
	List.WriteBinary(out);
}

/** \brief Write a Configuration in binary format to a .Config file with the given prefix. 
 *  \param prefix File name prefix (".Config" is appended).
 *  \param List   Configuration to write. */
inline void WriteBinary(std::string prefix, const Configuration & List)
{
	prefix+=std::string(".Config");
	std::fstream ofile(prefix.c_str(), std::fstream::out | std::fstream::binary);
	WriteBinary(ofile, List);
}

/** \brief Read a Configuration in binary format from an input stream. 
 *  \param ifile Input stream to read from.
 *  \return      The deserialized Configuration. */
inline Configuration ReadBinary(std::istream & ifile)
{
	return Configuration(ifile);
}

/** \brief Read a Configuration in binary format from a .Config file with the given prefix. 
 *  \param prefix File name prefix (".Config" is appended).
 *  \return       The deserialized Configuration. */
inline Configuration ReadBinary(std::string prefix)
{
	prefix+=std::string(".Config");
	std::fstream ifile(prefix.c_str(), std::fstream::in | std::fstream::binary);
	return ReadBinary(ifile);
}

/** \brief Read a Configuration from a .pos file with the given prefix. 
 *  \param prefix File name prefix (".pos" is appended).
 *  \return       The Configuration read from file. */
Configuration ReadPos(std::string prefix);

/** \brief Read a Configuration in .pos format from an input stream. 
 *  \param ifile Input stream in .pos format.
 *  \return      The Configuration read from the stream. */
Configuration ReadPos(std::istream & ifile);

/** \brief Read a Configuration from a HOOMD XML input stream.
 *  \param ifile Input stream in HOOMD XML format.
 *  \return      The Configuration read from the stream. */
Configuration ReadHoomdXml(std::istream & ifile);

/** \brief Read a Configuration from a HOOMD .xml file with the given prefix.
 *  \param prefix File name prefix (".xml" is appended).
 *  \return       The Configuration read from file. */
Configuration ReadHoomdXml(std::string prefix);

/** \brief Read a Configuration from a stealth-format input stream with the given side length.
 *  \param file       Input stream in stealth format.
 *  \param SideLength Side length used for the (cubic) box basis vectors.
 *  \return           The Configuration read from the stream. */
Configuration ReadStealthOutput(std::istream & file, double SideLength);

/** \brief Read a Configuration from a stealth-format input stream (unit side length).
 *  \param file Input stream in stealth format.
 *  \return     The Configuration read from the stream. */
Configuration ReadStealthOutput(std::istream & file);

/** \brief Read a Configuration from a GROMACS .gro input stream. 
 *  \param file Input stream in GROMACS .gro format.
 *  \return     The Configuration read from the stream. */
Configuration ReadGro(std::istream & file);

/** \brief Read a Configuration of relative coordinates from an input stream (unit side length).
 *  \param file Input stream of relative coordinates.
 *  \return     The Configuration read from the stream. */
Configuration ReadCoordinate(std::istream & file);

/** \brief Read a Configuration of relative coordinates with the given side length.
 *  \param file       Input stream of relative coordinates.
 *  \param SideLength Side length used for the (cubic) box basis vectors.
 *  \return           The Configuration read from the stream. */
Configuration ReadCoordinate(std::istream & file, double SideLength);

/** \brief Read a Configuration of Cartesian coordinates from an input stream (unit side length).
 *  \param file Input stream of Cartesian coordinates.
 *  \return     The Configuration read from the stream. */
Configuration ReadCartesianCoordinate(std::istream & file);

/** \brief Read a Configuration of Cartesian coordinates with the given side length.
 *  \param file       Input stream of Cartesian coordinates.
 *  \param SideLength Side length used for the (cubic) box basis vectors.
 *  \return           The Configuration read from the stream. */
Configuration ReadCartesianCoordinate(std::istream & file, double SideLength);

/** \brief Read a SpherePacking from a Torquato-Jiao (TJ) format input stream.
 *  \param inFile Input stream in TJ format.
 *  \return       The SpherePacking read from the stream. */
SpherePacking ReadTJOutput(std::istream & inFile);

/** \brief Read a SpherePacking from a TJ-format .dat file with the given prefix.
 *  \param prefix File name prefix (".dat" is appended).
 *  \return       The SpherePacking read from file. */
SpherePacking ReadTJOutput(const std::string & prefix);


/** \brief Write a SpherePacking in TJ format to an output file stream.
 *  \param ofile Output file stream to write to.
 *  \param pk    SpherePacking to write. */
void WriteTJOutput(std::fstream & ofile, const SpherePacking & pk);

/** \brief Write a SpherePacking in TJ format to a .dat file with the given prefix.
 *  \param prefix File name prefix (".dat" is appended).
 *  \param pk     SpherePacking to write. */
void WriteTJOutput(const std::string & prefix, const SpherePacking & pk);

/** \brief Output Cartesian coordinates of particles with the dimension and basis vectors
 *  \param conf       Configuration to write.
 *  \param OutputName File name prefix (".txt" is appended). */
//output Cartesian coordinates of particles with the dimension and basis vectors
void WriteConfiguration(const Configuration & conf, const std::string & OutputName);

/** \brief Write a SpherePacking (dimension, basis vectors, coordinates and radii) to a .txt file. 
 *  \param conf       SpherePacking to write.
 *  \param OutputName File name prefix (".txt" is appended). */
void WriteConfiguration(const SpherePacking & conf, const std::string & OutputName);

/** \brief Write a generic PeriodicCellList (dimension, basis vectors, coordinates) to a .txt file. 
 *  \param conf       PeriodicCellList to write.
 *  \param OutputName File name prefix (".txt" is appended). */
template<typename T> void WriteConfiguration(const PeriodicCellList<T> & conf, const std::string & OutputName);

/** \brief Read a Configuration from a .txt file written by WriteConfiguration. 
 *  \param conf      Configuration to populate.
 *  \param InputName File name prefix (".txt" is appended). */
void ReadConfiguration(Configuration & conf, const std::string & InputName);

/** \brief Read a SpherePacking from a .txt file written by WriteConfiguration. 
 *  \param conf      SpherePacking to populate.
 *  \param InputName File name prefix (".txt" is appended). */
void ReadConfiguration(SpherePacking & conf, const std::string & InputName);

#endif
