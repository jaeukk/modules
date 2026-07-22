/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	2022 (migrated to modules/ 2026-07) */

/** \file main.cpp
 * \brief Command-line driver that generates ensembles of selected disordered
 *        models and writes them as a ConfigurationPack or as per-configuration
 *        SpherePacking text files.
 *
 * Adapted from the legacy c++/Models/CLI.cpp against the modular `models` library.
 * The stdin token sequence is kept identical to the legacy CLI (see
 * c++/Models/README.md) so legacy generation scripts run unchanged:
 *   RSA:           dim N rho phi seed threads   nConfig prefix type
 *   HSF:           dim N rho phi seed threads interval   nConfig prefix type
 *   Poisson:       dim N rho mode               nConfig prefix type
 *   Perturbations: dim seed distribution [params]   nConfig inPrefix prefix type
 * Legacy tokens are mapped onto the module APIs via setters; where the module has
 * no equivalent (RSA seeds from a high-resolution clock) the token is consumed and
 * a warning is printed.
 *
 * Hardening over the legacy CLI: the model pointer is a nullptr-checked unique_ptr
 * (no indeterminate/leaked raw pointer); all interactive input is validated
 * (stream state, 1 <= d <= MaxDimension, positive bounded counts, finite positive
 * densities) and an invalid model aborts cleanly; SpherePacking output file names
 * are built by formatting only the index and concatenating the prefix as data (the
 * user prefix is never used as a printf format string); the unsupported
 * Perturbations + SpherePacking combination is rejected up front.
 */

#include <cstring>
#include <cstdio>
#include <string>
#include <sstream>
#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include <memory>
#include <algorithm>

/* in $(core) */
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>
#include <etc.h>

/* in $(models) */
#include <GenerateConfigs.h>
#include <RSA_Gen.h>
#include <HSF_Gen.h>
#include <GenPerturb.h>

/** \brief Library-wide verbosity level (defined by each executable per repo convention). */
size_t Verbosity = 2;

//------------------------------------------------------------------------------
// Validated interactive reads.  Each prints a prompt, echoes the value, and
// returns false (after a clear message) on a bad stream state or out-of-range
// value so the caller can abort the model cleanly.
//------------------------------------------------------------------------------

/** \brief Read and validate a spatial dimension in [1, MaxDimension].
 *  \param in   Input stream.  \param out  Output stream.  \param[out] d  Dimension.
 *  \return true on a valid read. */
static bool ReadDim(std::istream & in, std::ostream & out, int & d) {
	out << "Dimensions = ";
	if (!(in >> d)) { out << "invalid dimension input.\n"; return false; }
	out << d << "\n";
	if (d < 1 || d > (int)MaxDimension) {
		out << "dimension must be in [1, " << (int)MaxDimension << "].\n";
		return false;
	}
	return true;
}

/** \brief Read and validate a positive, reasonably-bounded integer count.
 *  \param in in-stream.  \param out out-stream.  \param label prompt label.  \param[out] n value.
 *  \return true on a valid read. */
static bool ReadCount(std::istream & in, std::ostream & out, const char * label, long long & n) {
	out << label << " = ";
	if (!(in >> n)) { out << "invalid " << label << ".\n"; return false; }
	out << n << "\n";
	if (n <= 0 || n > 100000000LL) {
		out << label << " must be a positive count (<= 1e8).\n";
		return false;
	}
	return true;
}

/** \brief Read and validate a finite positive real (e.g. number density).
 *  \param in in-stream.  \param out out-stream.  \param label prompt label.  \param[out] v value.
 *  \return true on a valid read. */
static bool ReadPosReal(std::istream & in, std::ostream & out, const char * label, double & v) {
	out << label << " = ";
	if (!(in >> v)) { out << "invalid " << label << ".\n"; return false; }
	out << v << "\n";
	if (!(v > 0.0) || !std::isfinite(v)) {
		out << label << " must be finite and positive.\n";
		return false;
	}
	return true;
}

/** \brief Read a plain integer (any value; used for seeds).
 *  \param in in-stream.  \param out out-stream.  \param label prompt label.  \param[out] v value.
 *  \return true on a valid read. */
static bool ReadInt(std::istream & in, std::ostream & out, const char * label, long long & v) {
	out << label << " = ";
	if (!(in >> v)) { out << "invalid " << label << ".\n"; return false; }
	out << v << "\n";
	return true;
}

/** \brief A fixed-N Poisson (or grand-canonical Binomial) point-pattern generator.
 *
 * The module PoissonPtGen fills a pre-sized box; for a from-scratch CLI we need a
 * generator that also builds its own cubic box of side (numPrts/rho)^{1/d}. This
 * thin ConfigGen restores the legacy "Poisson" CLI behavior against the module API.
 */
class PoissonCLI : public ConfigGen {
	int mode = 0;			//< 0 = canonical (fixed N), 1 = grand-canonical (Poisson-distributed N)
	RandomGenerator rng;	//< random number generator
public:
	/** \brief Construct from validated values (legacy Poisson tokens: dim N rho mode).
	 *  \param dimension      Spatial dimension.
	 *  \param numberOfPrts   Number of particles.
	 *  \param numberDensity  Number density (sets the cubic box side to (N/rho)^{1/d}).
	 *  \param canonicalMode  0 = canonical, 1 = grand-canonical. */
	PoissonCLI(int dimension, long long numberOfPrts, double numberDensity, int canonicalMode) : ConfigGen() {
		this->d = (DimensionType)dimension;
		this->numPrts = (size_t)numberOfPrts;
		this->rho = numberDensity;
		this->mode = canonicalMode;
		this->rng.seed((int)GetPreciseClock());
	}

	/** \brief Reset the number density. \param r  New number density. */
	virtual void SetRho(double r) { this->rho = r; }

	/** \brief Generate a Poisson sphere packing (radius-0 points) in a cubic box.
	 *  \param[out] p  Output point pattern as a SpherePacking of radius-0 spheres. */
	virtual void GenerateP(SpherePacking & p) {
		double L = std::pow((double)this->numPrts / this->rho, 1.0 / this->d);
		std::vector<GeometryVector> basis(this->d, GeometryVector(static_cast<int>(this->d)));
		for (DimensionType i = 0; i < this->d; i++)
			basis[i].x[i] = L;
		p = SpherePacking((DimensionType)this->d, &basis[0], L);
		size_t N = (this->mode == 1) ? this->rng.PoissonRNG((double)this->numPrts) : this->numPrts;
		GeometryVector rel(static_cast<int>(this->d));
		for (size_t i = 0; i < N; i++) {
			for (DimensionType j = 0; j < this->d; j++)
				rel.x[j] = this->rng.RandomDouble();
			p.Insert(0.0, rel);
		}
	}

	/** \brief Generate a Poisson point configuration.
	 *  \param[out] c  Output point configuration. */
	virtual void GenerateC(Configuration & c) {
		SpherePacking p; this->GenerateP(p);
		c = Configuration(p, "a");
	}
};

/** \brief Helper that draws input configurations from a ConfigurationPack.
 *  Used by transformation models (e.g. Perturbations) that require a source ensemble. */
class ReadConfigPack
{
public:
	double Rescale;						//< isotropic rescaling factor applied to each loaded configuration
	size_t StartIdx;					//< index of the first configuration to read
	ConfigurationPack p;				//< the opened configuration pack
	RandomGenerator rng;				//< random generator (used for non-sequential sampling)
	bool sequential_sampling = true;	//< true: read consecutively; false: sample uniformly at random

	/** \brief Open a ConfigurationPack whose prefix is read from ifile.
	 *  \param ifile     Input stream supplying the pack prefix.
	 *  \param ofile     Output stream for prompts.
	 *  \param Rescale   Isotropic rescaling factor applied to each loaded configuration.
	 *  \param StartIdx  Index of the first configuration to read. */
	ReadConfigPack(std::istream & ifile, std::ostream & ofile, double Rescale, size_t StartIdx = 0)
	{
		this->Rescale = Rescale;
		this->StartIdx = StartIdx;
		ofile << "Input Prefix:";
		std::string prefix;
		ifile >> prefix;
		p.Open(prefix);
	}

	/** \brief Return the i-th source configuration (rescaled).
	 *  \param i  Offset from StartIdx (sequential) / ignored (random sampling).
	 *  \return   The loaded and rescaled configuration. */
	Configuration operator() (size_t i)
	{
		Configuration result;
		if (sequential_sampling) {
			result = p.GetConfig(i + this->StartIdx);
		}
		else {
			result = p.GetConfig((int)(std::floor(p.NumConfig() * rng.RandomDouble())));
		}
		result.Rescale(Rescale);
		return result;
	}
};

/** \brief CLI entry point.
 *  \param argc  Argument count.
 *  \param argv  Argument vector; `-h` lists models, `models` runs the generator dialog.
 *  \return      0 on normal exit, 1 on an input/usage error. */
int main(int argc, char ** argv) {
	std::vector<std::string> models{ "Poisson", "RSA", "HSF", "Perturbations" };

	/* input configuration (only used by transformation models) */
	bool input_required = false;
	std::function<const Configuration(size_t i)> GetConfigsFunction = nullptr;
	double Rescale = 1.0;
	size_t StartIdx = 0;

	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	std::unique_ptr<ConfigGen> model;	//owns the selected model; nullptr until one is chosen

	if (argc <= 1) {
		ofile << "Usage: " << argv[0] << " -h | models\n";
		return 1;
	}

	if (strcmp(argv[1], "-h") == 0) {
		ofile << "Available models\n";
		for (auto iter = models.begin(); iter != models.end(); ++iter)
			ofile << *iter << "\t";
		ofile << "\n";
		return 0;
	}
	else if (strcmp(argv[1], "models") == 0) {
		std::string modelName;
		ofile << "Model name = ";
		if (!(ifile >> modelName)) { ofile << "No model name given.\n"; return 1; }

		if (std::find(models.begin(), models.end(), modelName) == models.end()) {
			ofile << modelName << " is not in the available model list.\n";
			return 1;
		}

		// ---- Construct the model, reading tokens in the legacy order ------------
		if (modelName == "Poisson") {
			// legacy tokens: dim N rho mode
			int dim; long long N; double rho; long long mode;
			if (!ReadDim(ifile, ofile, dim)) return 1;
			if (!ReadCount(ifile, ofile, "Particle number", N)) return 1;
			if (!ReadPosReal(ifile, ofile, "Number density", rho)) return 1;
			if (!ReadInt(ifile, ofile, "Mode (0 = canonical/ 1 = grand canonical)", mode)) return 1;
			if (mode != 0 && mode != 1) { ofile << "mode must be 0 or 1.\n"; return 1; }
			model.reset(new PoissonCLI(dim, N, rho, (int)mode));
		}
		else if (modelName == "RSA") {
			// legacy tokens: dim N rho phi seed threads
			int dim; long long N; double rho, phi; long long seed, threads;
			if (!ReadDim(ifile, ofile, dim)) return 1;
			if (!ReadCount(ifile, ofile, "Particle number", N)) return 1;
			if (!ReadPosReal(ifile, ofile, "Number density", rho)) return 1;
			if (!ReadPosReal(ifile, ofile, "Packing fraction", phi)) return 1;
			if (!ReadInt(ifile, ofile, "Random seed", seed)) return 1;
			if (!ReadCount(ifile, ofile, "Number of threads", threads)) return 1;
			std::unique_ptr<RSA_Gen> rsa(new RSA_Gen(dim, (int)N, phi, rho));
			rsa->SetNumThreads((int)threads);
			//RSA_Gen seeds from a high-resolution clock and has no seed setter; consume + warn.
			ofile << "[note] RSA_Gen ignores the seed token (" << seed
				  << "): it seeds from a high-resolution clock.\n";
			model = std::move(rsa);
		}
		else if (modelName == "HSF") {
			// legacy tokens: dim N rho phi seed threads interval
			int dim; long long N; double rho, phi; long long seed, threads, interval;
			if (!ReadDim(ifile, ofile, dim)) return 1;
			if (!ReadCount(ifile, ofile, "Particle number", N)) return 1;
			if (!ReadPosReal(ifile, ofile, "Number density", rho)) return 1;
			if (!ReadPosReal(ifile, ofile, "Packing fraction", phi)) return 1;
			if (!ReadInt(ifile, ofile, "Random seed", seed)) return 1;
			if (!ReadCount(ifile, ofile, "Number of threads", threads)) return 1;
			if (!ReadCount(ifile, ofile, "Sampling interval", interval)) return 1;
			std::unique_ptr<HSF_Gen> hsf(new HSF_Gen(dim, (int)N, phi, rho));
			hsf->setting.SetRandomSeed((int)seed);
			hsf->setting.SetNumThreads((int)threads);
			hsf->SetInterval((int)interval);
			model = std::move(hsf);
		}
		else if (modelName == "Perturbations") {
			// legacy tokens: dim seed distribution [params]
			int dim; long long seed;
			if (!ReadDim(ifile, ofile, dim)) return 1;
			if (!ReadInt(ifile, ofile, "Random seed", seed)) return 1;
			std::string dist;
			ofile << "Type of distribution = ";
			if (!(ifile >> dist)) { ofile << "invalid distribution.\n"; return 1; }
			ofile << dist << "\n";
			std::string dl = dist;
			std::transform(dl.begin(), dl.end(), dl.begin(), [](unsigned char ch){ return (char)std::tolower(ch); });
			int nparam = (dl == "uniform" || dl == "gaussian") ? 1
					   : ((dl == "powerlaw" || dl == "power-law") ? 2 : 0);
			//Re-assemble the tokens GenPerturb(stream) expects (dim, distribution, params)
			//so its own parser configures the distribution; the seed is applied via SetSeed.
			std::stringstream distStream;
			distStream << dim << ' ' << dist;
			for (int k = 0; k < nparam; k++) {
				double param;
				ofile << "Distribution parameter " << k << " = ";
				if (!(ifile >> param)) { ofile << "invalid distribution parameter.\n"; return 1; }
				ofile << param << "\n";
				distStream << ' ' << param;
			}
			std::ostringstream sink;	//suppress GenPerturb's own duplicate prompts
			std::unique_ptr<GenPerturb> gp(new GenPerturb(distStream, sink));
			gp->SetSeed((int)seed);
			model = std::move(gp);
			input_required = true;
		}

		if (model == nullptr) {
			ofile << "No model was constructed. Aborting.\n";
			return 1;
		}

		long long NumConfigSigned = 0;
		if (!ReadCount(ifile, ofile, "Number of configurations", NumConfigSigned)) return 1;
		size_t NumConfig = (size_t)NumConfigSigned;

		if (input_required) {
			ReadConfigPack c(ifile, ofile, Rescale, StartIdx);
			//A ConfigurationPack knows its own size, so clamp NumConfig to what is available.
			c.sequential_sampling = true;
			std::cout << "The config pack contains " << c.p.NumConfig() << " configurations.";
			if ((size_t)c.p.NumConfig() < StartIdx) {
				GetConfigsFunction = nullptr;
				NumConfig = 0;
				std::cout << "StartIdx exceeds the size of ConfigPack\n";
			}
			else {
				size_t num = c.p.NumConfig() - StartIdx;
				if (NumConfig > num) {
					std::cout << "Set NumConfig to the pack size because the request exceeds it.\n";
					NumConfig = num;
				}
			}
			GetConfigsFunction = c;
		}

		/* Output ... */
		std::string output_name, output_type;

		ofile << "output prefix = ";
		Echo(ifile, ofile, output_name);

		ofile << "output type (ConfigPack/SpherePacking[.txt]) = ";
		Echo(ifile, ofile, output_type);

		//Reject the unsupported combination up front: a transformation model reads its
		//source as point Configurations, so SpherePacking output would write empty packings.
		if (input_required && output_type == "SpherePacking") {
			ofile << "Perturbations produces point configurations; SpherePacking output is "
					 "not supported for this model. Use ConfigPack.\n";
			return 1;
		}

		if (output_type == "ConfigPack") {
			ConfigurationPack save(output_name);

			ofile << "output file = " << output_name << ".ConfigPack\n";
			ofile << "Generating Configurations ... \n";
			progress_display pd(NumConfig);
			if (input_required && GetConfigsFunction != nullptr) {
				for (size_t i = 0; i < NumConfig; i++) {
					Configuration c = GetConfigsFunction(i);
					model->GenerateC(c);
					save.AddConfig(c);
					pd++;
				}
			}
			else {
				for (size_t i = 0; i < NumConfig; i++) {
					Configuration c;
					model->GenerateC(c);
					save.AddConfig(c);
					pd++;
				}
			}
		}
		else if (output_type == "SpherePacking") {
			ofile << "output file = " << output_name << "__%05zu.txt\n";
			progress_display pd(NumConfig);

			for (size_t i = 0; i < NumConfig; i++) {
				//Format ONLY the index with a literal format, then concatenate the prefix
				//as data (never treat the user-supplied prefix as a printf format string).
				char idxbuf[32];
				snprintf(idxbuf, sizeof(idxbuf), "__%05zu", i);
				std::string fname = output_name + idxbuf;
				SpherePacking c;
				model->GenerateP(c);
				WriteConfiguration(c, fname);
				pd++;
			}
		}
		else {
			ofile << "Unrecognized output type \"" << output_type << "\".\n";
			return 1;
		}
		ofile << "\n Done!\n";
		return 0;
	}
	else {
		ofile << "Unrecognized command \"" << argv[1] << "\". Use -h for the model list.\n";
		return 1;
	}
}
