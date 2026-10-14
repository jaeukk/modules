#ifndef MSPB_PACK_H
#define MSPB_PACK_H

#include <SpherePacking.h>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

// A lossless MSPB v1 snapshot. IDs are metadata, not array indices.
struct MSPBParticle {
    std::uint64_t id = 0;
    std::uint32_t species = 0;
    GeometryVector relative;
    double radius = 0.0;
};

struct MSPBSnapshot {
    std::uint64_t speciesCount = 0; // Includes declared but empty species.
    std::vector<GeometryVector> basis;
    std::vector<MSPBParticle> particles;

    // Reject invalid geometry, nonfinite values, duplicate IDs and bad labels/radii.
    void Validate() const;
    void Write(std::ostream &out) const;
    // Reads one snapshot from a seekable stream; trailing bytes are an error.
    static MSPBSnapshot Read(std::istream &in);
    // Geometry/radii view for core algorithms. IDs and labels remain in this snapshot.
    SpherePacking ToSpherePacking() const;
};

// Append-only, single-writer archive: <prefix>.MSPBPack. MSPB payloads are unchanged.
// Open scans framing into an in-memory offset index; GetConfig validates its payload.
// Exceptions report invalid indices, malformed data, stale handles, and I/O failures.
// Not a checkpoint or a concurrent-writer store. Reopen after another handle writes.
class MSPB_Pack {
public:
    MSPB_Pack() = default;
    explicit MSPB_Pack(const std::string &prefix) { Open(prefix); }
    void Open(const std::string &prefix); // Missing file is an empty, not-yet-created pack.
    std::size_t NumConfig() const { return entries_.size(); }
    void AddConfig(const MSPBSnapshot &snapshot);
    MSPBSnapshot GetConfig(std::size_t index) const;
    void ImportFile(const std::string &filename);
    void ExportFile(std::size_t index, const std::string &filename) const;
    void Clear(); // Explicitly discard every record, retaining a valid empty archive.
    const std::string &FileName() const { return filename_; }
private:
    struct Entry { std::uint64_t offset, length; };
    std::string filename_;
    std::vector<Entry> entries_;
    std::uint64_t fileSize_ = 0;
    void CheckCurrent() const;
};

#endif
