#include "MSPB_Pack.h"
#include "MultiComponentPairStat.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <omp.h>

size_t Verbosity = 0;

namespace {
std::size_t count(const std::string &text) {
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
        throw std::invalid_argument("Expected a non-negative integer: " + text);
    const auto value = std::stoull(text);
    if (value > std::numeric_limits<std::size_t>::max()) throw std::out_of_range("Index too large");
    return static_cast<std::size_t>(value);
}
double number(const std::string &text) {
    std::size_t end;
    const double value = std::stod(text, &end);
    if (end != text.size()) throw std::invalid_argument("Invalid numeric parameter");
    return value;
}
std::string indexed(const std::string &prefix, std::size_t index, const char *suffix) {
    std::ostringstream out;
    out << prefix << '_' << std::setw(4) << std::setfill('0') << index << suffix;
    return out.str();
}
int usage() {
    std::cerr << "Usage:\n"
              << "  mspb_pack.out pack <pack prefix> <snapshot.mspb> [snapshot.mspb ...]\n"
              << "  mspb_pack.out inspect <pack prefix>\n"
              << "  mspb_pack.out unpack <pack prefix> <start> <count> <output prefix>\n"
              << "  mspb_pack.out pairstat <pack prefix> <start> <count> <report prefix> Kmax dK Rmax dR samples seed\n";
    return 1;
}
}

int main(int argc, char **argv) {
    try {
        if (argc < 3) return usage();
        const std::string command(argv[1]);
        if (command == "pack") {
            if (argc < 4) return usage();
            MSPB_Pack pack(argv[2]);
            for (int i = 3; i < argc; ++i) pack.ImportFile(argv[i]);
            std::cout << "# snapshots = " << pack.NumConfig() << '\n';
            return 0;
        }
        if (command != "inspect" && command != "unpack" && command != "pairstat") return usage();
        if ((command == "inspect" && argc != 3) || (command == "unpack" && argc != 6) ||
            (command == "pairstat" && argc != 12)) return usage();
        MSPB_Pack pack(argv[2]);
        if (!std::filesystem::exists(pack.FileName())) throw std::runtime_error("Archive does not exist");
        if (command == "inspect") {
            std::cout << "# snapshots = " << pack.NumConfig() << '\n';
            for (std::size_t i = 0; i < pack.NumConfig(); ++i) {
                const auto snapshot = pack.GetConfig(i);
                std::cout << i << ": dimension=" << snapshot.basis.size()
                          << " particles=" << snapshot.particles.size()
                          << " species=" << snapshot.speciesCount << '\n';
            }
            return 0;
        }
        const auto start = count(argv[3]), n = count(argv[4]);
        if (start > pack.NumConfig() || n > pack.NumConfig() - start)
            throw std::out_of_range("Requested range exceeds archive size");
        MultiComponentPairStatOptions options;
        if (command == "pairstat") {
            options.kMax = number(argv[6]); options.kBinWidth = number(argv[7]);
            options.rMax = number(argv[8]); options.rStep = number(argv[9]);
            options.windowSamples = count(argv[10]);
            std::size_t end;
            options.seed = std::stoi(argv[11], &end);
            if (end != std::string(argv[11]).size()) throw std::invalid_argument("Invalid seed");
            options.Validate();
            omp_set_num_threads(1); // Stable offline comparison; no optimizer state involved.
        }
        // Detect output collisions before writing any file in the selected range.
        const char *suffix = command == "unpack" ? ".mspb" : ".txt";
        for (std::size_t i = start; i < start + n; ++i)
            if (std::filesystem::exists(indexed(argv[5], i, suffix)))
                throw std::runtime_error("Output already exists: " + indexed(argv[5], i, suffix));
        for (std::size_t i = start; i < start + n; ++i) {
            const auto path = indexed(argv[5], i, suffix);
            if (command == "unpack") {
                pack.ExportFile(i, path);
            } else {
                const auto snapshot = pack.GetConfig(i);
                if (snapshot.speciesCount > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
                    throw std::invalid_argument("Species count exceeds PairStat adapter range");
                if (snapshot.basis.size() != 2 && snapshot.basis.size() != 3)
                    throw std::invalid_argument("PairStat adapter supports 2D and 3D");
                if (snapshot.particles.empty()) throw std::invalid_argument("Cannot analyze an empty packing");
                MultiStealthyPotential potential(static_cast<DimensionType>(snapshot.basis.size()), 0.0);
                potential.SetPacking(snapshot.ToSpherePacking());
                potential.SetNumSpecies(static_cast<std::size_t>(snapshot.speciesCount));
                for (std::size_t p = 0; p < snapshot.particles.size(); ++p)
                    potential.SetSpeciesOfParticle(p, snapshot.particles[p].species);
                const auto results = ComputeMultiComponentPairStatistics(potential, options);
                std::ofstream out(path);
                WriteMultiComponentPairStatistics(out, results, options);
                out.close();
                if (!out) throw std::runtime_error("Could not finish report: " + path);
            }
        }
        std::cout << n << " snapshots processed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
