#include "MSPB_Pack.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace {
static_assert(CHAR_BIT == 8 && sizeof(double) == 8 &&
              std::numeric_limits<double>::is_iec559, "MSPB requires binary64 doubles");
const std::string snapshotMagic("MSPB\r\n\x1a\n", 8);
const std::string packMagic("MSPP\r\n\x1a\n", 8);

void fail(const std::string &message) { throw std::runtime_error("MSPB: " + message); }
void readBytes(std::istream &in, char *data, std::size_t size) {
    if (!in.read(data, static_cast<std::streamsize>(size))) fail("truncated input");
}
std::uint64_t integer(std::istream &in, unsigned width) {
    unsigned char bytes[8] = {};
    readBytes(in, reinterpret_cast<char *>(bytes), width);
    std::uint64_t value = 0;
    for (unsigned b = 0; b < width; ++b) value |= std::uint64_t(bytes[b]) << (8*b);
    return value;
}
void integer(std::ostream &out, std::uint64_t value, unsigned width) {
    char bytes[8];
    for (unsigned b = 0; b < width; ++b) {
        const unsigned char octet = static_cast<unsigned char>(value & 255);
        std::memcpy(bytes + b, &octet, 1);
        value >>= 8;
    }
    out.write(bytes, width);
}
double real(std::istream &in) {
    const auto bits = integer(in, 8);
    double value;
    std::memcpy(&value, &bits, 8);
    return value;
}
void real(std::ostream &out, double value) {
    std::uint64_t bits;
    std::memcpy(&bits, &value, 8);
    integer(out, bits, 8);
}
void magic(std::istream &in, const std::string &expected) {
    char bytes[8];
    readBytes(in, bytes, 8);
    if (std::string(bytes, 8) != expected) fail("invalid signature");
}
std::uint64_t remaining(std::istream &in) {
    const auto pos = in.tellg();
    if (pos == std::streampos(-1)) fail("input must be seekable");
    in.seekg(0, std::ios::end);
    const auto end = in.tellg();
    if (end == std::streampos(-1) || end < pos) fail("cannot determine input size");
    in.seekg(pos);
    if (!in) fail("cannot seek input");
    return static_cast<std::uint64_t>(end - pos);
}
void packHeader(std::ostream &out) {
    out.write(packMagic.data(), 8);
    integer(out, 1, 4);
    integer(out, 0, 4);
}
}

void MSPBSnapshot::Validate() const {
    const auto dim = basis.size();
    if (dim == 0 || dim > MaxDimension) fail("unsupported dimension");
    // Pivoted elimination rejects singular boxes before core cell-list construction.
    long double matrix[MaxDimension][MaxDimension] = {};
    for (std::size_t d = 0; d < dim; ++d) {
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
        if (basis[d].Dimension != dim) fail("inconsistent basis dimension");
#endif
        for (std::size_t k = 0; k < dim; ++k) {
            if (!std::isfinite(basis[d].x[k])) fail("nonfinite basis");
            matrix[d][k] = basis[d].x[k];
        }
    }
    for (std::size_t k = 0; k < dim; ++k) {
        std::size_t pivot = k;
        for (std::size_t d = k + 1; d < dim; ++d)
            if (std::abs(matrix[d][k]) > std::abs(matrix[pivot][k])) pivot = d;
        if (matrix[pivot][k] == 0 || !std::isfinite(matrix[pivot][k])) fail("singular basis");
        for (std::size_t j = 0; j < dim; ++j) std::swap(matrix[k][j], matrix[pivot][j]);
        for (std::size_t d = k + 1; d < dim; ++d) {
            const auto factor = matrix[d][k] / matrix[k][k];
            for (std::size_t j = k + 1; j < dim; ++j) matrix[d][j] -= factor * matrix[k][j];
        }
    }
    std::unordered_set<std::uint64_t> ids;
    for (const auto &p : particles) {
        if (p.species >= speciesCount) fail("species label outside declared range");
        if (!ids.insert(p.id).second) fail("duplicate particle ID");
        if (!std::isfinite(p.radius) || p.radius <= 0) fail("invalid radius");
#ifdef GEOMETRYVECTOR_RECORDDIMENSION
        if (p.relative.Dimension != dim) fail("inconsistent coordinate dimension");
#endif
        for (std::size_t d = 0; d < dim; ++d)
            if (!std::isfinite(p.relative.x[d])) fail("nonfinite coordinate");
    }
}

void MSPBSnapshot::Write(std::ostream &out) const {
    Validate(); // No output for an invalid in-memory snapshot.
    out.write(snapshotMagic.data(), 8);
    integer(out, 1, 4);
    integer(out, basis.size(), 4);
    integer(out, particles.size(), 8);
    integer(out, speciesCount, 8);
    integer(out, 3, 4);
    for (const auto &b : basis)
        for (std::size_t d = 0; d < basis.size(); ++d) real(out, b.x[d]);
    for (const auto &p : particles) {
        integer(out, p.id, 8);
        integer(out, p.species, 4);
        for (std::size_t d = 0; d < basis.size(); ++d) real(out, p.relative.x[d]);
        real(out, p.radius);
    }
    if (!out) fail("output failed");
}

MSPBSnapshot MSPBSnapshot::Read(std::istream &in) {
    magic(in, snapshotMagic);
    if (integer(in, 4) != 1) fail("unsupported snapshot version");
    const auto dim = integer(in, 4);
    const auto count = integer(in, 8);
    MSPBSnapshot result;
    result.speciesCount = integer(in, 8);
    if (integer(in, 4) != 3) fail("unsupported snapshot flags");
    if (dim == 0 || dim > MaxDimension) fail("unsupported dimension");
    const auto bytes = remaining(in), basisBytes = 8*dim*dim, recordBytes = 20 + 8*dim;
    // Check against actual file length before allocating or multiplying a count.
    if (bytes < basisBytes || (bytes - basisBytes) % recordBytes != 0 ||
        count != (bytes - basisBytes) / recordBytes || count > result.particles.max_size())
        fail("particle count does not match file size");
    for (std::size_t d = 0; d < dim; ++d) {
        GeometryVector b(static_cast<int>(dim));
        for (std::size_t k = 0; k < dim; ++k) b.x[k] = real(in);
        result.basis.push_back(b);
    }
    result.particles.reserve(static_cast<std::size_t>(count));
    for (std::uint64_t i = 0; i < count; ++i) {
        MSPBParticle p;
        p.id = integer(in, 8);
        p.species = static_cast<std::uint32_t>(integer(in, 4));
        p.relative = GeometryVector(static_cast<int>(dim));
        for (std::size_t d = 0; d < dim; ++d) p.relative.x[d] = real(in);
        p.radius = real(in);
        result.particles.push_back(p);
    }
    result.Validate();
    return result;
}

SpherePacking MSPBSnapshot::ToSpherePacking() const {
    Validate();
    auto vectors = basis;
    // One cell per box direction avoids huge cell grids for large physical boxes.
    double cellSize = 0;
    for (const auto &b : basis) {
        double length = 0;
        for (std::size_t d = 0; d < basis.size(); ++d) length = std::hypot(length, b.x[d]);
        cellSize = std::max(cellSize, length);
    }
    if (!std::isfinite(cellSize) || cellSize <= 0) fail("box scale unsupported by core");
    SpherePacking result(static_cast<DimensionType>(basis.size()), vectors.data(), cellSize);
    if (!std::isfinite(result.PeriodicVolume()) || result.PeriodicVolume() <= 0)
        fail("box volume unsupported by core");
    for (const auto &p : particles) result.Insert(p.radius, p.relative);
    return result;
}

void MSPB_Pack::Open(const std::string &prefix) {
    if (prefix.empty()) throw std::invalid_argument("MSPB_Pack: empty prefix");
    const std::string path = prefix + ".MSPBPack";
    std::vector<Entry> entries;
    std::uint64_t size = 0;
    if (std::filesystem::exists(path)) {
        std::ifstream in(path, std::ios::binary);
        if (!in) fail("cannot open archive " + path);
        size = remaining(in);
        magic(in, packMagic);
        if (integer(in, 4) != 1 || integer(in, 4) != 0) fail("unsupported archive version/flags");
        std::uint64_t offset = 16;
        while (offset < size) {
            if (size - offset < 8) fail("incomplete archive frame");
            const auto length = integer(in, 8);
            offset += 8;
            if (length < 36 || length > size - offset) fail("invalid archive frame length");
            entries.push_back({offset, length});
            offset += length;
            in.seekg(static_cast<std::streamoff>(offset));
            if (!in) fail("cannot seek archive");
        }
    }
    filename_ = path;
    entries_ = std::move(entries);
    fileSize_ = size;
}

void MSPB_Pack::CheckCurrent() const {
    if (filename_.empty()) throw std::logic_error("MSPB_Pack: call Open first");
    const bool exists = std::filesystem::exists(filename_);
    if ((!exists && fileSize_ != 0) || (exists &&
        (fileSize_ == 0 || std::filesystem::file_size(filename_) != fileSize_)))
        fail("archive changed; reopen before use");
}

void MSPB_Pack::AddConfig(const MSPBSnapshot &snapshot) {
    std::ostringstream buffer(std::ios::out | std::ios::binary);
    snapshot.Write(buffer);
    const auto bytes = buffer.str();
    CheckCurrent();
    const std::uint64_t base = fileSize_ ? fileSize_ : 16;
    if (bytes.size() > static_cast<std::uint64_t>(std::numeric_limits<std::streamoff>::max()) - base - 8)
        fail("archive too large");
    // Reserve before writing so allocation failure cannot leave an unindexed record.
    if (entries_.size() == entries_.capacity()) {
        if (entries_.size() == entries_.max_size()) fail("too many archive records");
        entries_.reserve(entries_.size() + std::min(entries_.size() + 1,
                                                  entries_.max_size() - entries_.size()));
    }
    std::ofstream out(filename_, std::ios::binary | std::ios::app);
    if (!out) fail("cannot append archive " + filename_);
    if (!fileSize_) packHeader(out);
    integer(out, bytes.size(), 8);
    out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    out.close();
    if (!out) {
        std::error_code error;
        if (fileSize_) std::filesystem::resize_file(filename_, fileSize_, error);
        else std::filesystem::remove(filename_, error);
        fail(error ? "append failed; rollback failed: " + error.message() : "append failed; rolled back");
    }
    entries_.push_back({base + 8, bytes.size()});
    fileSize_ = base + 8 + bytes.size();
}

MSPBSnapshot MSPB_Pack::GetConfig(std::size_t index) const {
    CheckCurrent();
    if (index >= entries_.size()) throw std::out_of_range("MSPB_Pack: index out of range");
    const auto &entry = entries_[index];
    std::ifstream in(filename_, std::ios::binary);
    if (!in) fail("cannot read archive");
    in.seekg(static_cast<std::streamoff>(entry.offset));
    if (entry.length > std::string().max_size() ||
        entry.length > static_cast<std::uint64_t>(std::numeric_limits<std::streamsize>::max()))
        fail("snapshot too large");
    std::string bytes(static_cast<std::size_t>(entry.length), '\0');
    readBytes(in, &bytes[0], bytes.size());
    std::istringstream snapshot(bytes, std::ios::in | std::ios::binary);
    return MSPBSnapshot::Read(snapshot);
}

void MSPB_Pack::ImportFile(const std::string &filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) fail("cannot open snapshot " + filename);
    AddConfig(MSPBSnapshot::Read(in));
}

void MSPB_Pack::ExportFile(std::size_t index, const std::string &filename) const {
    const auto snapshot = GetConfig(index);
    // Do not accidentally truncate an archive or an existing result.
    if (std::filesystem::exists(filename)) fail("export destination already exists: " + filename);
    std::ofstream out(filename, std::ios::binary);
    if (!out) fail("cannot create snapshot " + filename);
    snapshot.Write(out);
    out.close();
    if (!out) fail("cannot finish snapshot " + filename);
}

void MSPB_Pack::Clear() {
    CheckCurrent();
    std::ofstream out(filename_, std::ios::binary | std::ios::trunc);
    if (!out) fail("cannot clear archive " + filename_);
    packHeader(out);
    out.close();
    if (!out) fail("cannot finish clearing archive");
    entries_.clear();
    fileSize_ = 16;
}
