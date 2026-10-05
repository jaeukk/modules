#include <PeriodicCellList.h>
#include <cassert>
#include <string>

// Binary storage must preserve doubles exactly, independently of text conversion.
static void Check(const SpherePacking & a, const SpherePacking & b)
{
    assert(a.GetDimension() == b.GetDimension());
    assert(a.NumParticle() == b.NumParticle());
    for (DimensionType d = 0; d < a.GetDimension(); ++d) {
        for (DimensionType k = 0; k < a.GetDimension(); ++k)
            assert(a.GetBasisVector(d).x[k] == b.GetBasisVector(d).x[k]);
        for (size_t i = 0; i < a.NumParticle(); ++i)
            assert(a.GetRelativeCoordinates(i).x[d] == b.GetRelativeCoordinates(i).x[d]);
    }
    for (size_t i = 0; i < a.NumParticle(); ++i)
        assert(a.GetCharacteristics(i) == b.GetCharacteristics(i));
    assert(a.GetMaxRadius() == b.GetMaxRadius());
}

int main(int argc, char ** argv)
{
    assert(argc == 2);
    const std::string prefix(argv[1]);
    GeometryVector basis[3] = {GeometryVector(3), GeometryVector(3), GeometryVector(3)};
    basis[0].x[0] = 2.123456789012345;
    basis[1].x[0] = 0.25; basis[1].x[1] = 1.5;
    basis[2].x[1] = 0.125; basis[2].x[2] = 2.0;
    SpherePacking source(3, basis, 1.0), empty(3, basis, 1.0);
    GeometryVector position(3);
    position.x[0] = 0.1234567890123456; position.x[1] = 0.4; position.x[2] = 0.6;
    source.Insert(0.1234567890123456, position);
    position.x[0] = 0.75;
    source.Insert(0.2345678901234567, position);
    SpherePackingPack pack(prefix);
    pack.AddConfig(source);
    SpherePackingPack reopened(prefix);
    assert(reopened.NumConfig() == 1);
    Check(source, reopened.GetConfig(0));
    reopened.AddConfig(empty);
    SpherePackingPack appended(prefix);
    assert(appended.NumConfig() == 2);
    Check(empty, appended.GetConfig(1));
    Check(source, appended.GetConfig(0));
    appended.Clear();
    SpherePackingPack cleared(prefix);
    assert(cleared.NumConfig() == 0);
    cleared.AddConfig(source);
    Check(source, SpherePackingPack(prefix).GetConfig(0));
    std::cout << "PASS: exact SpherePackingPack binary restoration, reopen, append, clear\n";
}
