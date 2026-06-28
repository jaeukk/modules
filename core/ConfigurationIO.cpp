/** \file ConfigurationIO.cpp
 *  \brief Implementations of the Configuration/SpherePacking file-format
 *         readers and writers declared in ConfigurationIO.h.  Split out of
 *         PeriodicCellList.cpp for readability. */
#include "ConfigurationIO.h"
#include <exception>

/** \brief Write a Configuration to a .pos text file with the given prefix. (CLAUDE)
 *  \param prefix File name prefix (".pos" is appended).
 *  \param List   Configuration to write. */
void Output(std::string prefix, const Configuration & List)
{
	prefix+=std::string(".pos");
	std::fstream out(prefix.c_str(), std::fstream::out);
	out.precision(17);

	Output(out, List);
}

/** \brief Write a Configuration in .pos text format to an output stream. (CLAUDE)
 *  \param out  Output stream to write to.
 *  \param List Configuration to write. */
void Output(std::ostream & out, const Configuration & List)
{
	const Configuration * list = & List;
	//determine Rc of output
	//GeometryVector temp(list->GetDimension());
	//for( ::DimensionType i=0; i<list->GetDimension(); i++)
	//{
	//	GeometryVector t=list->GetBasisVector(i);
	//	if(temp.Dot(t)>0)
	//		temp.AddFrom(t);
	//	else
	//		temp.MinusFrom(t);
	//}
	//double Rc=2*std::sqrt(temp.Modulus2())/std::pow(list->NumParticle(), 1.0/list->GetDimension());
	double Rc=std::pow(100*list->PeriodicVolume()/list->NumParticle(), 1.0/list->GetDimension());
	//////////////////////////////////////
	GeometryVector origin(list->GetDimension());
	size_t count=0;
	list->IterateThroughNeighbors(origin, Rc, [&count](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void {count++;} );
	out<<count<<'\n';
	//the second line is a comment
	out<<'\n';
	//.xyz part, used by visualization sofwares
	list->IterateThroughNeighbors(origin, Rc, [&out, &list](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void 
	{
		out<<list->GetCharacteristics(SourceAtom).name<<" \t"; 
		shift.OutputCoordinate(out, list->GetDimension());
		out<<'\n';
	} );


	out<<"\n\n****************************************\n\n";

	out<<"Primitive vectors\n";
	for(::DimensionType i=0; i<list->GetDimension(); i++)
	{
		out<<"a(";
		out<<i;
		out<<") = \t";
		list->GetBasisVector(i).OutputCoordinate(out, list->GetDimension());
		out<<'\n';
	}

	out<<"\n\nVolume =  "<<list->PeriodicVolume()<<"\n\n";
	out<<"Reciprocal vectors\n";
	for(::DimensionType i=0; i<list->GetDimension(); i++)
	{
		out<<"b(";
		out<<i;
		out<<") = \t";
		list->GetReciprocalBasisVector(i).OutputCoordinate(out, list->GetDimension());
		out<<'\n';
	}
	out<<"\n\nBasis Vectors:\nAtom    Lattice Coordinates                Cartesian Coordinates\n\n";
	for(size_t i=0; i<list->NumParticle(); i++)
	{
		//Configuration::particle * now=list->GetParticle(i);
		out<<list->GetCharacteristics(i).name<<" \t";
		list->GetRelativeCoordinates(i).OutputCoordinate(out, list->GetDimension());
		list->GetCartesianCoordinates(i).OutputCoordinate(out, list->GetDimension());
		out<<'\n';
	}
}
/** \brief Read a Configuration from a .pos file with the given prefix. (CLAUDE)
 *  \param prefix File name prefix (".pos" is appended).
 *  \return       The Configuration read from file. */
Configuration ReadPos(std::string prefix)
{
	std::string name=prefix+std::string(".pos");
	std::fstream file(name.c_str(), std::fstream::in);
	return ReadPos(file);
}
/** \brief Read a Configuration in .pos format from an input stream. (CLAUDE)
 *  \param ifile Input stream in .pos format.
 *  \return      The Configuration read from the stream. */
Configuration ReadPos(std::istream & ifile)
{
	const std::string flag1("Primitive");
	const std::string flag2("Basis");

	std::string current;

	size_t NumParticle=0;

	::DimensionType dimension=0;
	do
	{
		std::getline(ifile, current);
		if(ifile.eof() || ifile.fail())
		{
			return Configuration(0, nullptr, 0.0, false);
		}
	}
	while(current.find(flag1)!=0);
	for(;;)
	{
		std::getline(ifile, current);
		if(current.size()<3)
			break;
		if(current[2]>='0' && current[2]<='9')
			dimension++;
		else
			break;
	}

	do
	{
		std::getline(ifile, current);
	}
	while(current.find(flag2)!=0);
	std::getline(ifile, current);
	std::getline(ifile, current);
	while(ifile.eof()==false)
	{
		char junk[1000];
		ifile.getline(junk, 1000);

		if(ifile.eof()==false)
			NumParticle++;
	}
	ifile.clear();
	ifile.seekg(0, std::ios::beg);


	do
	{
		std::getline(ifile, current);
	}
	while(current.find(flag1)!=0);

	std::vector<GeometryVector> base;
	for(int i=0; i<dimension; i++)
	{
		char a;
		do
		{
			a=ifile.get();
		}
		while(a!='=');
		base.push_back(GeometryVector(dimension));
		base.back().InputCoordinate(ifile, dimension);
	}
	double volume= ::Volume(&base[0], dimension);
	double cellsize = std::pow(volume/(NumParticle/ParticlePerCell), static_cast<double>(1)/dimension);
	Configuration result(dimension, &base[0], cellsize);

	do
	{
		std::getline(ifile, current);
	}
	while(current.find(flag2)!=0);
	std::getline(ifile, current);
	std::getline(ifile, current);
	while(ifile.eof()==false)
	{
		char name[128];
		ifile >> name;

		GeometryVector loc(dimension);
		loc.InputCoordinate(ifile, dimension);

		GeometryVector junk(dimension);
		junk.InputCoordinate(ifile, dimension);

		if(ifile.eof()==false)
			result.Insert(name, loc);
	}
	return result;
}
/** \brief Read a Configuration from a HOOMD .xml file with the given prefix. (CLAUDE)
 *  \param prefix File name prefix (".xml" is appended).
 *  \return       The Configuration read from file. */
Configuration ReadHoomdXml(std::string prefix)
{
	prefix+=std::string(".xml");
	std::fstream ifile(prefix.c_str());
	return ReadHoomdXml(ifile);
}
/** \brief Read a Configuration from a HOOMD XML input stream. (CLAUDE)
 *  \param ifile Input stream in HOOMD XML format.
 *  \return      The Configuration read from the stream. */
Configuration ReadHoomdXml(std::istream & ifile)
{
	std::string temp;
	std::string flag1("<configuration");
	std::string flag2("<type");

	DimensionType dim;
	size_t NumParticle;
	double volume=1;

	std::vector<GeometryVector> basis;
	std::vector<double> coords;

	do
	{
		std::getline(ifile, temp);
	}
	while(temp.find(flag1)!=0);

	{
		std::stringstream stream;
		stream<<temp;

		//read a number
		for(;;)
		{
			char c=stream.peek();
			if(c>='0' && c<='9')
				break;
			if(c=='.')
				break;
			if(stream.eof())
				break;
			stream>>c;
		}
		stream>>NumParticle;//this is actually timestep. junk it
		//read a number
		for(;;)
		{
			char c=stream.peek();
			if(c>='0' && c<='9')
				break;
			if(c=='.')
				break;
			if(stream.eof())
				break;
			stream>>c;
		}
		stream>>dim;
		//read a number
		for(;;)
		{
			char c=stream.peek();
			if(c>='0' && c<='9')
				break;
			if(c=='.')
				break;
			if(stream.eof())
				break;
			stream>>c;
		}
		stream>>NumParticle;
	}

	std::getline(ifile, temp);
	{
		std::stringstream stream;
		stream<<temp;
		double length;

		for(DimensionType d=0; d<dim; d++)
		{
			//read a number
			for(;;)
			{
				char c=stream.peek();
				if(c>='0' && c<='9')
					break;
				if(c=='.')
					break;
				if(stream.eof())
					break;
				stream>>c;
			}
			stream>>length;
			GeometryVector tempVec(dim);
			tempVec.x[d]=length;
			basis.push_back(tempVec);
			volume=volume*length;
		}
	}

	double cellsize = std::pow(volume/(NumParticle/ParticlePerCell), static_cast<double>(1)/dim);
	Configuration result(dim, &basis[0], cellsize);
	std::getline(ifile, temp);
	for(size_t i=0; i<dim*NumParticle; i++)
	{
		double tempD;
		ifile>>tempD;
		coords.push_back(tempD);
	}


	do
	{
		std::getline(ifile, temp);
	}
	while(temp.find(flag2)!=0);
	for(size_t i=0; i<NumParticle; i++)
	{
		std::getline(ifile, temp);
		assert(temp.length()<3);
		GeometryVector rel(dim);
		for(DimensionType j=0; j<dim; j++)
			rel.x[j]=coords[i*dim+j]/basis[j].x[j];

		result.Insert(temp.c_str(), rel);
	}

	return result;
}
/** \brief Read a Configuration of relative coordinates from an input stream (unit side length). (CLAUDE)
 *  \param file Input stream of relative coordinates.
 *  \return     The Configuration read from the stream. */
Configuration ReadCoordinate(std::istream & file)
{
	return ReadCoordinate(file, 1.0);
}
/** \brief Read a Configuration of relative coordinates with the given side length. (CLAUDE)
 *  \param file       Input stream of relative coordinates.
 *  \param SideLength Side length used for the (cubic) box basis vectors.
 *  \return           The Configuration read from the stream. */
Configuration ReadCoordinate(std::istream & file, double SideLength)
{
	DimensionType dim=0;
	size_t NumP=0;
	std::string temp;

	//determine dimension
	std::getline(file, temp);
	std::stringstream str;
	str<<temp;
	for(;;)
	{
		double dtemp=-1.0;
		str>>dtemp;
		if(dtemp!=-1.0)
			dim++;
		else
			break;
	}

	//determine NumP
	while(file.eof()==false)
	{
		std::getline(file, temp);
		NumP++;
	}
	file.clear();
	file.seekg(0, std::ios::beg);

	std::vector<GeometryVector> basis;
	for(DimensionType i=0; i<dim; i++)
	{
		basis.push_back(GeometryVector(dim));
		basis.back().x[i]=SideLength;
	}
	double CellSize=SideLength/std::pow(static_cast<double>(NumP), 1.0/static_cast<double>(dim));
	Configuration result(dim, &basis[0], CellSize);

	for(size_t i=0; i<NumP; i++)
	{
		GeometryVector rel(dim);
		for(size_t j=0; j<dim; j++)
			file>>rel.x[j];
		result.Insert("A", rel);
	}

	return result;
}
/** \brief Read a Configuration of Cartesian coordinates from an input stream (unit side length). (CLAUDE)
 *  \param file Input stream of Cartesian coordinates.
 *  \return     The Configuration read from the stream. */
Configuration ReadCartesianCoordinate(std::istream & file)
{
	return ReadCartesianCoordinate(file, 1.0);
}
/** \brief Read a Configuration of Cartesian coordinates with the given side length. (CLAUDE)
 *  \param file       Input stream of Cartesian coordinates.
 *  \param SideLength Side length used for the (cubic) box basis vectors.
 *  \return           The Configuration read from the stream. */
Configuration ReadCartesianCoordinate(std::istream & file, double SideLength)
{
	DimensionType dim=0;
	size_t NumP=0;
	std::string temp;

	//determine dimension
	std::getline(file, temp);
	std::stringstream str;
	str<<temp;
	for(;;)
	{
		double dtemp=-1.0;
		str>>dtemp;
		if(dtemp!=-1.0)
			dim++;
		else
			break;
	}

	//determine NumP
	while(file.eof()==false)
	{
		std::getline(file, temp);
		NumP++;
	}
	file.clear();
	file.seekg(0, std::ios::beg);

	std::vector<GeometryVector> basis;
	for(DimensionType i=0; i<dim; i++)
	{
		basis.push_back(GeometryVector(dim));
		basis.back().x[i]=SideLength;
	}
	double CellSize=SideLength/std::pow(static_cast<double>(NumP), 1.0/static_cast<double>(dim));
	Configuration result(dim, &basis[0], CellSize);

	for(size_t i=0; i<NumP; i++)
	{
		GeometryVector rel(dim);
		for(size_t j=0; j<dim; j++)
			file>>rel.x[j];
		rel=result.CartesianCoord2RelativeCoord(rel);
		result.Insert("A", rel);
	}

	return result;
}

/** \brief Read a Configuration from a GROMACS .gro input stream (reads only atoms named 'O'). (CLAUDE)
 *  \param file Input stream in GROMACS .gro format.
 *  \return     The Configuration read from the stream. */
Configuration ReadGro(std::istream & file)
{
	GeometryVector basis[ ::MaxDimension];
	DimensionType dim=3;
	size_t NumP=0;
	std::string temp;

	//determine NumP
	std::getline(file, temp);
	file>>NumP;
	for(size_t i=0; i<NumP; i++)
	{
		std::getline(file, temp);
	}
	std::getline(file, temp);
	//determine box size
	double Volume=1.0;
	for(int i=0; i<3; i++)
	{
		double l;
		file>>l;
		basis[i]=GeometryVector(3);
		basis[i].x[i]=l;
		Volume*=l;
	}
	file.clear();
	file.seekg(0, std::ios::beg);
	std::getline(file, temp);
	std::getline(file, temp);

	double CellSize=Volume/std::pow(static_cast<double>(NumP), 1.0/static_cast<double>(dim));
	Configuration result(dim, &basis[0], CellSize);

	for(size_t i=0; i<NumP; i++)
	{
		//Only read atoms named 'O'!
		file>>temp;
		file>>temp;
		if(temp.at(0)=='O')
		{
			file>>temp;
			GeometryVector rel(dim);
			for(size_t j=0; j<dim; j++)
				file>>rel.x[j];
			rel=result.CartesianCoord2RelativeCoord(rel);
			result.Insert("A", rel);
		}
		else
			std::getline(file, temp);
	}

	return result;
}
/** \brief Read a Configuration from a stealth-format input stream with the given side length. (CLAUDE)
 *  \param file       Input stream in stealth format.
 *  \param SideLength Side length used for the (cubic) box basis vectors.
 *  \return           The Configuration read from the stream. */
Configuration ReadStealthOutput(std::istream & file, double SideLength)
{
	DimensionType dim=0;
	size_t NumP=0;
	std::string temp;
	std::vector<double> OneOverSize;//the side length of the input file

	//determine dimension
	std::getline(file, temp);
	{
		std::stringstream str;
		str<<temp;
		str>>temp;
		str>>dim;
	}

	//determine size
	std::getline(file, temp);
	{
		std::stringstream str;
		str<<temp;
		str>>temp;
		for(DimensionType i=0; i<dim; i++)
		{
			double tt;
			str>>tt;
			OneOverSize.push_back(1.0/tt);
		}
	}

	std::getline(file, temp);
	//determine NumP
	while(file.eof()==false)
	{
		std::getline(file, temp);
		NumP++;
	}
	file.clear();
	file.seekg(0, std::ios::beg);
	std::getline(file, temp);
	std::getline(file, temp);

	std::vector<GeometryVector> basis;
	for(DimensionType i=0; i<dim; i++)
	{
		basis.push_back(GeometryVector(dim));
		basis.back().x[i]=SideLength;
	}
	double CellSize=SideLength/std::pow(static_cast<double>(NumP), 1.0/static_cast<double>(dim));
	Configuration result(dim, &basis[0], CellSize);

	for(size_t i=0; i<NumP; i++)
	{
		GeometryVector rel(dim);
		for(size_t j=0; j<dim; j++)
		{
			file>>rel.x[j];
			rel.x[j]*=OneOverSize[j];
		}
		result.Insert("A", rel);
	}

	return result;
}
/** \brief Read a Configuration from a stealth-format input stream (unit side length). (CLAUDE)
 *  \param file Input stream in stealth format.
 *  \return     The Configuration read from the stream. */
Configuration ReadStealthOutput(std::istream & file)
{
	return ReadStealthOutput(file, 1.0);
}

/** \brief Read a SpherePacking from a Torquato-Jiao (TJ) format input stream. (CLAUDE)
 *  \param inFile Input stream in TJ format.
 *  \return       The SpherePacking read from the stream. */
SpherePacking ReadTJOutput(std::istream & inFile)
{
	std::string input;
	std::getline(inFile,input);
	std::stringstream ss1(input);
	int dim, N;
	ss1 >> dim;
	std::string buffer;
	ss1 >> buffer; //"HS"
	std::string packType;
	ss1 >> packType; //"poly" or "mono"
	getline(inFile,input); //Line 2 is junk.
	inFile >> N;

	double radMono;
	if (packType == "mono") 
	{//This line is the diameter value for monodisperse packings
		inFile >> radMono;
		radMono /= 2.0; //Convert from diameter to radius
	}
	else if (packType != "poly")
	{
		std::cerr << "Error in ReadTJOutput : Unsupported packType!\n";
		return SpherePacking();
	}

	// lattice vectors next
	GeometryVector basis[ ::MaxDimension];
	for (int i=0; i<dim; i++) 
	{
		basis[i]=GeometryVector(dim);
		for (int j=0; j<dim; j++) 
		{
			double val;
			inFile >> val;
			basis[i].x[j]=val;
		}
	}
	SpherePacking packing(dim, basis, ::MaxDistance);

	// then more trash
	getline(inFile,input);
	getline(inFile,input);

	// then coordinates and radii
	for (int i=0; i<N; i++) 
	{
		double radii;
		GeometryVector car(dim);
		for (int d=0; d<dim; d++) 
		{
			inFile >> car.x[d];
		}
		if (packType == "poly") 
		{
			inFile >> radii;
		}
		else 
		{
			radii = radMono;
		}
		packing.Insert(radii, packing.CartesianCoord2RelativeCoord(car));
	}
	packing.SetCellSize(packing.GetMaxRadius()*1.05);

	return packing; //Got a packing!

}
/** \brief Read a SpherePacking from a TJ-format .dat file with the given prefix. (CLAUDE)
 *  \param prefix File name prefix (".dat" is appended).
 *  \return       The SpherePacking read from file. */
SpherePacking ReadTJOutput(const std::string & prefix)
{
	std::string name = prefix + std::string(".dat");
	std::fstream ifile(name.c_str(), std::fstream::in);
	return ReadTJOutput(ifile);
}
/** \brief Write a SpherePacking in TJ format to an output file stream. (CLAUDE)
 *  \param ofile Output file stream to write to.
 *  \param pk    SpherePacking to write. */
void WriteTJOutput(std::fstream & ofile, const SpherePacking & pk)
{
	ofile.precision(17);
	DimensionType dim = pk.GetDimension();
	size_t num = pk.NumParticle();
	if (dim == 0 || num == 0)
	{
		std::cerr << "Error in WriteTJOutput : invalid packing!\n";
		return;
	}
	//determine mono or poly
	bool mono = true;
	for(size_t i=1; i<num; i++)
		if (pk.GetCharacteristics(i) != pk.GetCharacteristics(0))
		{
			mono = false;
			break;
		}

	ofile << dim << "	hs	";
	if (mono)
		ofile << "mono";
	else
		ofile << "poly";
	ofile << "\n" << num << " 1\n" << num << '\n';
	if (mono)
		ofile << 2 * pk.GetCharacteristics(0) << '\n';

	for (DimensionType i = 0; i < dim; i++)
	{
		for (DimensionType j = 0; j < dim; j++)
			ofile << pk.GetBasisVector(i).x[j] << " \t";
		ofile << '\n';
	}
	ofile << "T T T";
	for (size_t n = 0; n < num; n++)
	{
		ofile << "\n";
		for (DimensionType j = 0; j < dim; j++)
			ofile << pk.GetCartesianCoordinates(n).x[j] << " \t";
		if (!mono)
			ofile << pk.GetCharacteristics(n);
	}
}
/** \brief Write a SpherePacking in TJ format to a .dat file with the given prefix. (CLAUDE)
 *  \param prefix File name prefix (".dat" is appended).
 *  \param pk     SpherePacking to write. */
void WriteTJOutput(const std::string & prefix, const SpherePacking & pk)
{
	std::string name = prefix + std::string(".dat");
	std::fstream ifile(name.c_str(), std::fstream::out);
	return WriteTJOutput(ifile, pk);
}



/** \brief Write a Configuration (dimension, basis vectors, Cartesian coordinates) to a .txt file. (CLAUDE)
 *  \param conf       Configuration to write.
 *  \param OutputName File name prefix (".txt" is appended). */
void WriteConfiguration(const Configuration & conf, const std::string & OutputName)
{
	// dimension 0	0	0
	// basis vectors ...
	// Cartessian coordinates
	std::string name = OutputName + std::string(".txt");
	std::fstream output(name.c_str(), std::fstream::out);
	
	output << std::setprecision(17);
	output << conf.GetDimension() << "\n";
	for (int i = 0; i < conf.GetDimension(); i++) {
		output <<conf.GetBasisVector(i);
	}
	for (size_t i = 0; i < conf.NumParticle(); i++) {
		output << conf.GetCartesianCoordinates(i);
	}
	output.close();
}

/** \brief Write a SpherePacking (dimension, basis vectors, coordinates and radii) to a .txt file. (CLAUDE)
 *  \param conf       SpherePacking to write.
 *  \param OutputName File name prefix (".txt" is appended). */
void WriteConfiguration(const SpherePacking & conf, const std::string & OutputName)
{
	//dimension	0	0	0
	//basisvector1
	// ...
	//Relative coordinates1	radius 1
	//
	//....

	std::string name = OutputName + std::string(".txt");
	std::fstream ofile(name.c_str(), std::fstream::out);
	ofile << std::setprecision(10);
	ofile << conf.GetDimension() << "\n";
	for (int i = 0; i < conf.GetDimension(); i++) {
		ofile << conf.GetBasisVector(i);
	}
	//coordinates
	for (size_t i = 0; i < conf.NumParticle(); i++) {
		GeometryVector temp(conf.GetCartesianCoordinates(i));
		ofile << temp.x[0]<<"\t"<<temp.x[1]<<"\t"<<temp.x[2]<<"\t"<<temp.x[3]<<"\t"<<conf.GetCharacteristics(i)<<"\n";
	}
	ofile.close();
}

/** \brief Write a generic PeriodicCellList (dimension, basis vectors, coordinates) to a .txt file. (CLAUDE)
 *  \param conf       PeriodicCellList to write.
 *  \param OutputName File name prefix (".txt" is appended). */
template<typename T> void WriteConfiguration(const PeriodicCellList<T> & conf, const std::string & OutputName)
{

	std::string name = OutputName + std::string(".txt");
	std::fstream ofile(name.c_str(), std::fstream::out);
	ofile << std::setprecision(10);
	ofile << conf.GetDimension() << "\n";
	for (int i = 0; i < conf.GetDimension(); i++) {
		ofile << conf.GetBasisVector(i);
	}
	//coordinates
	for (size_t i = 0; i < conf.NumParticle(); i++) {
		GeometryVector temp(conf.GetCartesianCoordinates(i));
		ofile << temp.x[0] << "\t" << temp.x[1] << "\t" << temp.x[2] << "\t" << temp.x[3] << "\t" << conf.GetCharacteristics(i) << "\n";
	}
	ofile.close();
}
/** \brief Read a Configuration from a .txt file written by WriteConfiguration (assume # columns = 4). (CLAUDE)
 *  \param conf      Configuration to populate.
 *  \param InputName File name prefix (".txt" is appended). */
//assume # columns = 4
void ReadConfiguration(Configuration & conf, const std::string & InputName) {
	
	std::string name = InputName + std::string(".txt");
	std::fstream ifile(name.c_str(), std::fstream::in);
	
	if (ifile.good() == false)
	{
		std::cerr << "Error in ReadFunction : input stream is not good!\n" << InputName << "\n";
		return;
	}

	
	DimensionType d = 0;
	GeometryVector * basis;
	// set dimension
	ifile >> d;
	basis = new GeometryVector[d];

	// set basis vectors
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			ifile >> basis[i].x[j];
		}
		ifile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		basis[i].SetDimension(d);
	}
	//initialize conf
	conf = Configuration(d, basis, 1.0);
	while (ifile.eof() == false)
	{
		char dump[255] = {};		
		//input particle positions
		GeometryVector temp(static_cast<int> (d));
		for (size_t i = 0; i < ::MaxDimension; i++) {
			if (i < d) {
				ifile >> temp.x[i];
			}
			else {
				ifile.getline(dump, 255); i = ::MaxDimension;
			}
		}

		if (ifile.eof() == false)
			conf.Insert("a", conf.CartesianCoord2RelativeCoord(temp));
	}

	//std::cout << "the number of particles = " << conf.NumParticle() << "\n";

}

/** \brief Read a SpherePacking from a .txt file written by WriteConfiguration (assume # of columns is 5). (CLAUDE)
 *  \param conf      SpherePacking to populate.
 *  \param InputName File name prefix (".txt" is appended). */
//assume # of columns is 5
void ReadConfiguration(SpherePacking & conf, const std::string & InputName) {
	std::string name = InputName + std::string(".txt");
	std::fstream ifile(name.c_str(), std::fstream::in);

	if (ifile.good() == false)
	{
		std::cerr << "Error in ReadFunction : input stream is not good!\n"<<InputName<<"\n";
		return;
	}
	
	DimensionType d = 0;
	GeometryVector * basis;
	// set dimension
	ifile >> d;
	basis = new GeometryVector[d];

	// set basis vectors
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			ifile >> basis[i].x[j];
		}
		ifile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		basis[i].SetDimension(d);
	}
	//initialize conf
	conf = SpherePacking(d, basis, 1);
	while (ifile.eof() == false)
	{
		char dump[255] = {};
		//input particle positions
		GeometryVector temp(static_cast<int> (d));
		double radius = 0;
		for (size_t i = 0; i < ::MaxDimension; i++) {
			if (i < d) {
				ifile >> temp.x[i];
			}
			else {
				ifile >> dump;
			}
		}
		ifile >> radius;

		if (ifile.eof() == false)
			conf.Insert(radius, conf.CartesianCoord2RelativeCoord(temp));
	}
	
	
	/*std::vector<GeometryVector> temp;
	ReadFunction(temp, InputName);

	DimensionType d = (DimensionType)temp[0].x[0];
	GeometryVector *basis;
	basis = new GeometryVector[d];
	for (int i = 0; i < d; i++) {
		basis[i] = temp[i + 1];
		basis[i].SetDimension(d);
	}

	conf = SpherePacking(d, basis, 1.0);
	int numPrt = (temp.size() - (1 + d)) / 2;
	for (int i = 0; i < numPrt; i++) {
		conf.Insert(temp[2 * i + d + 2].x[0], temp[2 * i + d + 1]);
	}

	std::cout << "the number of particles = " << conf.NumParticle() << "\n";*/

}


