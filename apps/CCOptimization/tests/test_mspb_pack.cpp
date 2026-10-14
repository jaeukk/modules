#include "MSPB_Pack.h"
#include "MultiStealthyPotential.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

void check(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename F> void mustThrow(const std::string &message, F operation) {
    try {
        operation();
    } catch (const std::exception &) {
        return;
    }
    throw std::runtime_error("Expected exception: " + message);
}

void integer(std::string &bytes, std::uint64_t value, unsigned width) {
    for (unsigned i = 0; i < width; ++i) {
        bytes.push_back(static_cast<char>(value & 0xff));
        value >>= 8;
    }
}

void real(std::string &bytes, double value) {
    std::uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    integer(bytes, bits, 8);
}

void replaceInteger(std::string &bytes, std::size_t offset, std::uint64_t value, unsigned width) {
    check(offset + width <= bytes.size(), "test fixture offset");
    for (unsigned i = 0; i < width; ++i) {
        bytes[offset + i] = static_cast<char>(value & 0xff);
        value >>= 8;
    }
}

std::string rawMSPB(const MSPBSnapshot &snapshot) {
    const std::size_t dim = snapshot.basis.size();
    std::string bytes("MSPB\r\n\x1a\n", 8);
    integer(bytes, 1, 4);
    integer(bytes, dim, 4);
    integer(bytes, snapshot.particles.size(), 8);
    integer(bytes, snapshot.speciesCount, 8);
    integer(bytes, 3, 4);
    for (const auto &basis : snapshot.basis)
        for (std::size_t axis = 0; axis < dim; ++axis) real(bytes, basis.x[axis]);
    for (const auto &particle : snapshot.particles) {
        integer(bytes, particle.id, 8);
        integer(bytes, particle.species, 4);
        for (std::size_t axis = 0; axis < dim; ++axis) real(bytes, particle.relative.x[axis]);
        real(bytes, particle.radius);
    }
    return bytes;
}

std::string archive(const std::string &payload) {
    std::string bytes("MSPP\r\n\x1a\n", 8);
    integer(bytes, 1, 4);
    integer(bytes, 0, 4);
    integer(bytes, payload.size(), 8);
    bytes += payload;
    return bytes;
}

void writeFile(const std::string &name, const std::string &bytes) {
    std::ofstream out(name, std::ios::binary | std::ios::trunc);
    check(bool(out), "open test file for writing: " + name);
    out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    out.close();
    check(bool(out), "write test file: " + name);
}

std::string readFile(const std::string &name) {
    std::ifstream in(name, std::ios::binary);
    check(bool(in), "open test file for reading: " + name);
    std::ostringstream out;
    out << in.rdbuf();
    check(!in.bad(), "read test file: " + name);
    return out.str();
}

MSPBSnapshot sample(unsigned dim) {
    MSPBSnapshot snapshot;
    snapshot.speciesCount = 4; // Species 1 and 3 deliberately have no particles.
    for (unsigned vector = 0; vector < dim; ++vector) {
        GeometryVector basis(static_cast<int>(dim));
        basis.x[vector] = 3.0 + vector;
        if (vector == 1) basis.x[0] = 0.625;
        if (vector == 2) basis.x[1] = -0.25;
        snapshot.basis.push_back(basis);
    }
    for (unsigned i = 0; i < 3; ++i) {
        MSPBParticle particle;
        particle.id = 700 + 17 * i; // Nonsequential IDs are metadata.
        particle.species = i == 1 ? 2 : 0;
        particle.relative = GeometryVector(static_cast<int>(dim));
        for (unsigned axis = 0; axis < dim; ++axis)
            particle.relative.x[axis] = (i + axis + 1) / 16.0;
        particle.radius = (i + 1) / 8.0;
        snapshot.particles.push_back(particle);
    }
    return snapshot;
}

void compare(const MSPBSnapshot &actual, const MSPBSnapshot &expected) {
    check(actual.speciesCount == expected.speciesCount, "species count changed");
    check(actual.basis.size() == expected.basis.size(), "dimension changed");
    check(actual.particles.size() == expected.particles.size(), "particle count changed");
    for (std::size_t i = 0; i < expected.basis.size(); ++i)
        for (std::size_t axis = 0; axis < expected.basis.size(); ++axis)
            check(actual.basis[i].x[axis] == expected.basis[i].x[axis], "basis changed");
    for (std::size_t i = 0; i < expected.particles.size(); ++i) {
        const auto &a = actual.particles[i];
        const auto &e = expected.particles[i];
        check(a.id == e.id && a.species == e.species && a.radius == e.radius,
              "particle metadata changed");
        for (std::size_t axis = 0; axis < expected.basis.size(); ++axis)
            check(a.relative.x[axis] == e.relative.x[axis], "relative coordinate changed");
    }
}

MSPBSnapshot parse(const std::string &bytes) {
    std::istringstream in(bytes, std::ios::binary);
    return MSPBSnapshot::Read(in);
}

void invalidSnapshotTests(const MSPBSnapshot &valid) {
    auto invalid = valid;
    invalid.particles[1].id = invalid.particles[0].id;
    mustThrow("duplicate ID", [&] { invalid.Validate(); });
    invalid = valid;
    invalid.particles[0].species = 4;
    mustThrow("out-of-range species", [&] { invalid.Validate(); });
    invalid = valid;
    invalid.particles[0].radius = 0.0;
    mustThrow("zero radius", [&] { invalid.Validate(); });
    invalid.particles[0].radius = -1.0;
    mustThrow("negative radius", [&] { invalid.Validate(); });
    invalid.particles[0].radius = std::numeric_limits<double>::infinity();
    mustThrow("nonfinite radius", [&] { invalid.Validate(); });
    invalid = valid;
    invalid.particles[0].relative.x[0] = std::numeric_limits<double>::quiet_NaN();
    mustThrow("nonfinite relative coordinate", [&] { invalid.Validate(); });
    invalid = valid;
    invalid.basis[0].x[0] = std::numeric_limits<double>::infinity();
    mustThrow("nonfinite basis", [&] { invalid.Validate(); });
    invalid = valid;
    invalid.basis[1] = invalid.basis[0];
    mustThrow("singular basis", [&] { invalid.Validate(); });
    invalid = valid;
    invalid.speciesCount = 0;
    mustThrow("zero species with particles", [&] { invalid.Validate(); });

    std::ostringstream output;
    invalid = valid;
    invalid.particles[0].radius = 0.0;
    mustThrow("writer rejects invalid snapshot", [&] { invalid.Write(output); });
    check(output.str().empty(), "invalid writer modified stream");
}

void invalidPayloadTests(const MSPBSnapshot &valid) {
    const std::string good = rawMSPB(valid);
    auto bad = good;
    bad[0] = 'X';
    mustThrow("MSPB magic", [&] { parse(bad); });
    bad = good;
    replaceInteger(bad, 8, 2, 4);
    mustThrow("MSPB version", [&] { parse(bad); });
    bad = good;
    replaceInteger(bad, 12, 0, 4);
    mustThrow("zero dimension", [&] { parse(bad); });
    bad = good;
    replaceInteger(bad, 12, MaxDimension + 1, 4);
    mustThrow("excessive dimension", [&] { parse(bad); });
    bad = good;
    replaceInteger(bad, 16, std::numeric_limits<std::uint64_t>::max(), 8);
    mustThrow("impossible particle count", [&] { parse(bad); });
    bad = good;
    replaceInteger(bad, 24, 0, 8);
    mustThrow("invalid species count", [&] { parse(bad); });
    bad = good;
    replaceInteger(bad, 32, 7, 4);
    mustThrow("unknown MSPB flags", [&] { parse(bad); });
    bad = good;
    replaceInteger(bad, 32, 1, 4);
    mustThrow("missing required MSPB flag", [&] { parse(bad); });
    bad = good;
    bad.pop_back();
    mustThrow("truncated MSPB", [&] { parse(bad); });
    bad = good + "x";
    mustThrow("trailing MSPB byte", [&] { parse(bad); });
    bad = good;
    // First particle label follows the 36-byte header, D*D basis doubles, and ID.
    replaceInteger(bad, 36 + 8 * valid.basis.size() * valid.basis.size() + 8, 4, 4);
    mustThrow("invalid serialized label", [&] { parse(bad); });
    bad = good;
    const std::size_t firstRadius = 36 + 8 * valid.basis.size() * valid.basis.size()
                                    + 8 + 4 + 8 * valid.basis.size();
    replaceInteger(bad, firstRadius, 0, 8);
    mustThrow("invalid serialized radius", [&] { parse(bad); });
}

void snapshotTests(unsigned dim) {
    const auto expected = sample(dim);
    expected.Validate();
    const auto bytes = rawMSPB(expected);
    check(bytes.size() == 36 + 8 * dim * dim + 3 * (20 + 8 * dim), "MSPB fixture size");
    compare(parse(bytes), expected);
    std::ostringstream out(std::ios::binary);
    expected.Write(out);
    check(out.str() == bytes, "snapshot writer differs from independent MSPB fixture");
    const auto packing = expected.ToSpherePacking();
    check(packing.GetDimension() == dim && packing.NumParticle() == expected.particles.size(),
          "packing shape changed");
    for (std::size_t i = 0; i < expected.particles.size(); ++i) {
        check(packing.GetCharacteristics(i) == expected.particles[i].radius, "packing radius changed");
        for (unsigned axis = 0; axis < dim; ++axis)
            check(packing.GetRelativeCoordinates(i).x[axis] == expected.particles[i].relative.x[axis],
                  "packing relative coordinate changed");
    }
    invalidSnapshotTests(expected);
    invalidPayloadTests(expected);
}

void archiveTests() {
    const auto first = sample(2);
    const auto second = sample(3);
    const auto firstBytes = rawMSPB(first);
    const auto secondBytes = rawMSPB(second);
    MSPB_Pack pack("test_archive");
    check(pack.FileName() == "test_archive.MSPBPack", "archive filename");
    check(pack.NumConfig() == 0, "missing archive should be empty");
    std::ifstream absent(pack.FileName(), std::ios::binary);
    check(!absent, "opening missing archive created a file");
    mustThrow("empty archive index", [&] { pack.GetConfig(0); });
    pack.AddConfig(first);
    check(pack.NumConfig() == 1, "append count");
    check(readFile(pack.FileName()) == archive(firstBytes), "archive framing or payload changed");
    compare(pack.GetConfig(0), first);
    pack.ImportFile("second.mspb");
    check(pack.NumConfig() == 2, "import count");
    compare(pack.GetConfig(1), second);
    auto expected = archive(firstBytes);
    integer(expected, secondBytes.size(), 8);
    expected += secondBytes;
    check(readFile(pack.FileName()) == expected, "import did not preserve raw payload");
    pack.ExportFile(0, "exported.mspb");
    check(readFile("exported.mspb") == firstBytes, "export changed MSPB bytes");
    pack.ExportFile(1, "exported3.mspb");
    check(readFile("exported3.mspb") == secondBytes, "3D export changed MSPB bytes");
    mustThrow("out-of-range get", [&] { pack.GetConfig(2); });
    mustThrow("out-of-range export", [&] { pack.ExportFile(2, "unused.mspb"); });

    auto invalid = first;
    invalid.particles[0].radius = -1.0;
    mustThrow("invalid append", [&] { pack.AddConfig(invalid); });
    check(pack.NumConfig() == 2 && readFile(pack.FileName()) == expected,
          "invalid append altered archive");
    writeFile("bad-import.mspb", firstBytes + "x");
    mustThrow("invalid import", [&] { pack.ImportFile("bad-import.mspb"); });
    check(pack.NumConfig() == 2 && readFile(pack.FileName()) == expected,
          "invalid import altered archive");

    MSPB_Pack reopened("test_archive");
    check(reopened.NumConfig() == 2, "reopen count");
    compare(reopened.GetConfig(1), second);
    pack.AddConfig(first);
    mustThrow("stale reader", [&] { reopened.GetConfig(0); });
    mustThrow("stale writer", [&] { reopened.AddConfig(second); });
    reopened.Open("test_archive");
    check(reopened.NumConfig() == 3, "reopen after append");
    reopened.Clear();
    check(reopened.NumConfig() == 0, "clear count");
    check(readFile(reopened.FileName()) == std::string("MSPP\r\n\x1a\n\x01\0\0\0\0\0\0\0", 16),
          "clear did not retain valid empty header");
    reopened.AddConfig(second);
    check(reopened.NumConfig() == 1, "append after clear");
    compare(MSPB_Pack("test_archive").GetConfig(0), second);
}

void exporterCompatibilityTests() {
    for (unsigned dim : {2u, 3u}) {
        auto expected = sample(dim);
        GeometryVector basis[3] = {GeometryVector(static_cast<int>(dim)),
                                   GeometryVector(static_cast<int>(dim)),
                                   GeometryVector(static_cast<int>(dim))};
        for (unsigned i = 0; i < dim; ++i) basis[i] = expected.basis[i];
        Configuration configuration(static_cast<DimensionType>(dim), basis, 1.0);
        std::vector<double> radii;
        for (std::size_t i = 0; i < expected.particles.size(); ++i) {
            configuration.Insert("a", expected.particles[i].relative);
            expected.particles[i].id = i; // The existing exporter assigns global indices.
            radii.push_back(expected.particles[i].radius);
        }
        MultiStealthyPotential potential(static_cast<DimensionType>(dim), 1.0);
        potential.SetConfiguration(configuration);
        potential.SetNumSpecies(expected.speciesCount);
        for (std::size_t i = 0; i < expected.particles.size(); ++i)
            potential.SetSpeciesOfParticle(i, expected.particles[i].species);
        potential.SetSecureRadii(radii);
        std::ostringstream out(std::ios::binary);
        potential.ExportConfiguration(out, MultiStealthyPotential::ExportFormat::Binary);
        const auto bytes = out.str();
        check(bytes == rawMSPB(expected), "existing exporter differs from MSPB v1 fixture");
        compare(parse(bytes), expected);
        const auto file = "exporter-" + std::to_string(dim) + ".mspb";
        const auto prefix = "exporter-" + std::to_string(dim);
        writeFile(file, bytes);
        MSPB_Pack pack(prefix);
        pack.ImportFile(file);
        compare(pack.GetConfig(0), expected);
        pack.ExportFile(0, "roundtrip-" + file);
        check(readFile("roundtrip-" + file) == bytes,
              "archive changed existing exporter payload");
    }
}

void malformedArchiveTests() {
    const std::string good = archive(rawMSPB(sample(2)));
    auto expectBadOpen = [&](const std::string &name, const std::string &bytes) {
        writeFile(name + ".MSPBPack", bytes);
        mustThrow(name, [&] { MSPB_Pack pack(name); });
    };
    auto bad = good;
    bad[0] = 'X';
    expectBadOpen("bad-magic", bad);
    bad = good;
    replaceInteger(bad, 8, 2, 4);
    expectBadOpen("bad-version", bad);
    bad = good;
    replaceInteger(bad, 12, 1, 4);
    expectBadOpen("bad-flags", bad);
    expectBadOpen("short-header", good.substr(0, 15));
    expectBadOpen("short-length", good.substr(0, 20));
    bad = good;
    replaceInteger(bad, 16, std::numeric_limits<std::uint64_t>::max(), 8);
    expectBadOpen("huge-length", bad);
    expectBadOpen("short-payload", good.substr(0, good.size() - 1));
    expectBadOpen("trailing-byte", good + "x");

    // Framing is valid, so Open indexes the record; GetConfig checks MSPB contents.
    bad = rawMSPB(sample(2));
    bad[0] = 'X';
    writeFile("bad-payload.MSPBPack", archive(bad));
    MSPB_Pack pack("bad-payload");
    check(pack.NumConfig() == 1, "Open should index framed payload");
    mustThrow("invalid framed payload", [&] { pack.GetConfig(0); });
}

} // namespace

int main() {
    snapshotTests(2);
    snapshotTests(3);
    writeFile("second.mspb", rawMSPB(sample(3)));
    archiveTests();
    exporterCompatibilityTests();
    malformedArchiveTests();
    std::cout << "MSPB snapshot and pack tests passed.\n";
}
