#include "MultiStealthyPotential.h"
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <cstdint>
#include <cstring>

using Format = MultiStealthyPotential::ExportFormat;

static void check(bool ok) { if (!ok) throw std::runtime_error("export test failed"); }
template<class F> static void mustThrow(F f) {
    bool threw = false;
    try { f(); } catch (const std::exception &) { threw = true; }
    check(threw);
}

int main() {
    for (DimensionType dim : {2, 3}) {
      for (bool skewed : {false, true}) {
        GeometryVector basis[3] = {GeometryVector(dim), GeometryVector(dim), GeometryVector(dim)};
        for (DimensionType d = 0; d < dim; ++d) basis[d].x[d] = 4.0 + d;
        if (skewed) basis[1].x[0] = 1.0 / 7.0;
        Configuration c(dim, basis, 1.0);
        for (size_t i = 0; i < 3; ++i) {
            GeometryVector s(dim);
            for (DimensionType d = 0; d < dim; ++d) s.x[d] = (i + d + 1) / 11.0;
            c.Insert("a", s);
        }
        MultiStealthyPotential p(dim, 1.0);
        std::ostringstream out;
        mustThrow([&] { p.ExportConfiguration(out); });
        mustThrow([&] { p.ExportConfiguration(out, Format::Binary); });
        check(out.str().empty());
        p.SetConfiguration(c);
        p.SetNumSpecies(2);
        mustThrow([&] { p.ExportConfiguration(out); });
        mustThrow([&] { p.ExportConfiguration(out, Format::Binary); });
        check(out.str().empty());
        // Deliberately not i % species; exercise reassignment too.
        p.SetSpeciesOfParticle(0, 0);
        p.SetSpeciesOfParticle(0, 1);
        p.SetSpeciesOfParticle(1, 1);
        p.SetSpeciesOfParticle(2, 0);
        p.SetSecureRadii({0.125, 0.25, 0.375});
        out << std::scientific << std::setprecision(3);
        const auto flags = out.flags();
        const auto precision = out.precision();
        p.ExportConfiguration(out);
        check(out.flags() == flags && out.precision() == precision);
        std::istringstream input(out.str());
        std::string header;
        std::getline(input, header);
        check(header == "# dimension " + std::to_string(static_cast<unsigned int>(dim)));
        std::getline(input, header);
        check(header == "# particles 3");
        GeometryVector restoredBasis[3] = {GeometryVector(dim), GeometryVector(dim), GeometryVector(dim)};
        for (DimensionType d = 0; d < dim; ++d) {
            std::getline(input, header);
            std::istringstream line(header);
            std::string comment, label;
            line >> comment >> label;
            check(comment == "#" && label == "basis_vector-" + std::to_string(static_cast<unsigned int>(d)));
            for (DimensionType axis = 0; axis < dim; ++axis) {
                line >> restoredBasis[d].x[axis];
                check(bool(line) && restoredBasis[d].x[axis] == c.GetBasisVector(d).x[axis]);
            }
            check(!(line >> label));
        }
        Configuration restored(dim, restoredBasis, 1.0);
        std::getline(input, header);
        check(header == "# coordinates relative");
        std::getline(input, header);
        check(header == "# radius_units cartesian");
        std::getline(input, header);
        check(header == (dim == 2 ? "# id\tspecies\tsx\tsy\tradius" : "# id\tspecies\tsx\tsy\tsz\tradius"));
        for (size_t i = 0; i < 3; ++i) {
            size_t id; int species; double value;
            input >> id >> species;
            check(id == i && species == (i < 2 ? 1 : 0));
            GeometryVector relative(dim);
            for (DimensionType d = 0; d < dim; ++d) {
                input >> value;
                check(value == c.GetRelativeCoordinates(i).x[d]);
                relative.x[d] = value;
            }
            restored.Insert("a", relative);
            for (DimensionType d = 0; d < dim; ++d)
                check(restored.GetCartesianCoordinates(i).x[d] == c.GetCartesianCoordinates(i).x[d]);
            input >> value;
            check(value == (i + 1) * 0.125);
        }
        check(bool(input));
        std::string extra;
        check(!(input >> extra));
        std::ostringstream explicitText;
        p.ExportConfiguration(explicitText, Format::Text);
        check(explicitText.str() == out.str());

        std::ostringstream binary;
        binary << std::scientific << std::setprecision(3);
        const auto binaryFlags = binary.flags();
        p.ExportConfiguration(binary, Format::Binary);
        check(binary.flags() == binaryFlags && binary.precision() == 3);
        const std::string bytes = binary.str();
        check(bytes.size() == 36 + 8 * dim * dim + 3 * (20 + 8 * dim));
        check(bytes.substr(0, 8) == std::string("MSPB\r\n\x1a\n", 8));
        check(bytes.substr(8, 4) == std::string("\x01\0\0\0", 4));
        size_t offset = 8;
        auto readInteger = [&](size_t width) {
            check(offset + width <= bytes.size());
            std::uint64_t value = 0;
            for (size_t b = 0; b < width; ++b)
                value |= std::uint64_t(static_cast<unsigned char>(bytes[offset++])) << (8 * b);
            return value;
        };
        auto readReal = [&]() {
            const std::uint64_t bits = readInteger(8);
            double value;
            std::memcpy(&value, &bits, sizeof(value));
            return value;
        };
        check(readInteger(4) == 1);
        check(readInteger(4) == dim);
        check(readInteger(8) == 3);
        check(readInteger(8) == 2);
        check(readInteger(4) == 3);
        // Known little-endian binary64 encoding of the first basis component, 4.0.
        check(bytes.substr(offset, 8) == std::string("\0\0\0\0\0\0\x10\x40", 8));
        for (DimensionType d = 0; d < dim; ++d)
            for (DimensionType axis = 0; axis < dim; ++axis) {
                restoredBasis[d].x[axis] = readReal();
                check(restoredBasis[d].x[axis] == c.GetBasisVector(d).x[axis]);
            }
        Configuration binaryRestored(dim, restoredBasis, 1.0);
        for (size_t i = 0; i < 3; ++i) {
            check(readInteger(8) == i);
            check(readInteger(4) == (i < 2 ? 1 : 0));
            GeometryVector relative(dim);
            for (DimensionType d = 0; d < dim; ++d) {
                relative.x[d] = readReal();
                check(relative.x[d] == c.GetRelativeCoordinates(i).x[d]);
            }
            check(readReal() == (i + 1) * 0.125);
            binaryRestored.Insert("a", relative);
            for (DimensionType d = 0; d < dim; ++d)
                check(binaryRestored.GetCartesianCoordinates(i).x[d] == c.GetCartesianCoordinates(i).x[d]);
        }
        check(offset == bytes.size());
        std::ostringstream repeated;
        p.ExportConfiguration(repeated, Format::Binary);
        check(repeated.str() == bytes);
        std::ostringstream invalid;
        mustThrow([&] { p.ExportConfiguration(invalid, static_cast<Format>(99)); });
        check(invalid.str().empty());
        std::ostringstream failed;
        failed.setstate(std::ios::badbit);
        mustThrow([&] { p.ExportConfiguration(failed); });
        mustThrow([&] { p.ExportConfiguration(failed, Format::Binary); });
        // Export must reflect a refreshed final position, not the initial snapshot.
        GeometryVector moved(dim);
        moved.x[0] = 0.75;
        c.MoveParticle(0, moved);
        p.SetConfiguration(c);
        std::ostringstream refreshed;
        p.ExportConfiguration(refreshed);
        check(refreshed.str().find("0\t1\t0.75\t") != std::string::npos);
      }
    }
    std::cout << "Text and MSPB exporter tests passed (rectangular/skewed, 2D and 3D).\n";
}
