#include "Voronoi.h"
#include <DisplaySpheres.h>
#include <PairCorrelation.h>
// (legacy included etc1.h here, but Voronoi.cpp uses no etc1 symbol — dropped.)
#ifdef USE_MATHGL
#include <mgl2/mgl.h>
void PlotDisk2D(mglGraph & graph, double x, double y, double r, const char * MGLOptions, size_t NumLineSegments = 200);
#endif

//#ifdef USE_Voro++
#ifdef USE_CGAL

#include <CGAL/Cartesian_d.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Kernel_d/Sphere_d.h>
#include <CGAL/Delaunay_d.h>

//#define USE_CGAL_PRECISE

typedef CGAL::Cartesian_d<CGAL::Gmpq> Kernele;
#ifdef USE_CGAL_PRECISE 
typedef CGAL::Homogeneous_d<CGAL::Gmpq> Kerneli;
#else
typedef CGAL::Homogeneous_d<double> Kerneli;
#endif
typedef CGAL::Sphere_d<Kerneli> Sphere;
typedef CGAL::Convex_hull_d<Kernele> ConvexHull;
typedef ConvexHull::Point_d cPoint;
typedef ConvexHull::Vertex_handle cVertex_handle;
typedef ConvexHull::Simplex_handle cSimplex_handle;
typedef CGAL::Delaunay_d<Kerneli> Delaunay;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Simplex_handle Simplex_handle;
typedef Delaunay::Point_d Point;

template <typename point> GeometryVector CGALPoint2GeometryVector(const point & p, DimensionType d)
{
	GeometryVector result(d);
	for (DimensionType i = 0; i<d; i++)
		result.x[i] = CGAL::to_double(p.cartesian(i));
	return result;
}

// test if point "tested" is on the convex hull
bool IsOnConvexHull(const GeometryVector & tested, const std::vector<GeometryVector> & neighbors, DimensionType d)
{
	ConvexHull ch(d);
	cPoint p(d, tested.x, tested.x + d);
	cVertex_handle hTested = ch.insert(p);
	for (auto iter = neighbors.begin(); iter != neighbors.end(); iter++)
		ch.insert(cPoint(d, iter->x, iter->x + d));
	for (ConvexHull::Hull_vertex_iterator iter = ch.hull_vertices_begin(); iter != ch.hull_vertices_end(); iter++)
		if (iter->point() == hTested->point())
			return true;
	return false;
}

//calculate all vertexes of the voronoi cell of particle i, in vertices
//if pLinkage is not nullptr, it will be filled with linkage information of the format:
// (*pLinkage)[i*(Dimension+1)+j] is the index of a vertex linked with vertex i
// (*pLinkage)[i*(Dimension+1)+j] is -1 if vertex i have less than j linked vertexes.
void GetVoronoiCell(const PeriodicCellList<Empty> & c, size_t i, std::vector<GeometryVector> & vertices, std::vector<signed long> * pLinkage, std::set<size_t> * pNeighbors)
{
	DimensionType d = c.GetDimension();
	double radius = std::pow(c.PeriodicVolume() / c.NumParticle(), 1.0 / d)*1.8;
	double FurthestVertexDist = ::MaxDistance;
	GeometryVector ci = c.GetCartesianCoordinates(i);
	GeometryVector ri = c.GetRelativeCoordinates(i);
	struct VoronoiVertexInfo
	{
		Simplex_handle hDelaunaySimplex;
		size_t VerticesVectorLocation;
	};
	std::vector<VoronoiVertexInfo> vInfo;
	vertices.clear();

	//debug temp
	//size_t NumIteration = 0;
	while (2 * FurthestVertexDist>radius || vertices.empty())
	{
		//debug temp
		//NumIteration++;

		vInfo.clear();
		vertices.clear();
		radius *= 1.25;
		if (pNeighbors != nullptr)
			pNeighbors->clear();
		//Get neighbors within this radius
		std::vector<GeometryVector> neighbors;
		c.IterateThroughNeighbors(ri, radius, [&i, &ci, &neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
		{
			if (Sourceparticle == i && shift.Modulus2() == 0)
				return;
			neighbors.push_back(ci + shift);
		});
		//if no neighbors are found, then obviously radius is too small
		if (neighbors.empty())
			continue;
		//convex hull test, if the center particle is on convex hull then we are sure radius is not enough.
		if (IsOnConvexHull(ci, neighbors, d))
			continue;
		//Get Delaunay of neighbors
		Delaunay De(d);
		Point p(d, ci.x, ci.x + d);
		Vertex_handle hvi = De.insert(p);
		for (auto iter = neighbors.begin(); iter != neighbors.end(); iter++)
		{
			Point p(d, iter->x, iter->x + d);
			De.insert(p);
		}
		//for each simplex in Delaunay, if one of its vertex is hvi, then calculate the circumsphere center.
		std::list<Simplex_handle> simplices = De.all_simplices();
		for (auto iter = simplices.begin(); iter != simplices.end(); iter++)
		{
			bool has_hvi = false;
			bool IsValid = true;
			std::vector<Point> points;
			for (size_t i = 0; i <= d; i++)
			{
				Vertex_handle hv = De.vertex_of_simplex(*iter, i);
				if (hv == hvi)
					has_hvi = true;
				if (hv != nullptr)
					points.push_back(De.associated_point(hv));
				else
					IsValid = false;
			}
			if (IsValid == false)
				continue;
			if (has_hvi)
			{
				//calculate sphere center
				Sphere s(d, points.begin(), points.end());
				if (s.is_legal())
				{
					VoronoiVertexInfo temp;
					temp.hDelaunaySimplex = *iter;
					temp.VerticesVectorLocation = vertices.size();
					vInfo.push_back(temp);
					Point p_center = s.center();
					GeometryVector g_center = CGALPoint2GeometryVector(p_center, d);
					vertices.push_back(g_center);
				}

				if (pNeighbors != nullptr)
					for (int ii = 0; ii <= d; ii++)
					{
						GeometryVector car = CGALPoint2GeometryVector(De.associated_point(De.vertex_of_simplex(*iter, ii)), d);
						GeometryVector rel = c.CartesianCoord2RelativeCoord(car);
						bool found = false;
						size_t n = c.NumParticle();
						c.IterateThroughNeighbors(rel, LengthPrecision*0.01, [&found, &n](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
						{
							n = Sourceparticle;
							found = true;
						}, &found);

						//debug temp
						//std::cout << "found neighbor " << n << '\n';

						if (n < c.NumParticle())
						{
							if (n != i)
								pNeighbors->insert(n);
						}
						else
							std::cerr << "Error in GetVoronoiCell : predicted particle not found!\n";
					}
			}
		}
		//calculate FurthestVertexDist
		FurthestVertexDist = 0;
		for (auto iter = vertices.begin(); iter != vertices.end(); iter++)
		{
			GeometryVector dist = *iter - ci;
			if (FurthestVertexDist<dist.Modulus2())
				FurthestVertexDist = dist.Modulus2();
		}
		FurthestVertexDist = std::sqrt(FurthestVertexDist);

		if (vertices.empty() == false && 2 * FurthestVertexDist<radius && pLinkage != nullptr)
		{
			//calculate linkage
			pLinkage->resize((d + 1)*vertices.size(), -1);
			std::sort(vInfo.begin(), vInfo.end(), [](const VoronoiVertexInfo & left, const VoronoiVertexInfo & right) ->bool
			{
				return left.hDelaunaySimplex<right.hDelaunaySimplex;
			});
			for (auto iter = vInfo.begin(); iter != vInfo.end(); iter++)
			{
				size_t NeighborCount = 0;
				for (DimensionType i = 0; i <= d; i++)
				{
					Simplex_handle hNeighbor = De.opposite_simplex(iter->hDelaunaySimplex, i);
					auto NewFaceIterator = std::lower_bound(vInfo.begin(), vInfo.end(), hNeighbor, [](const VoronoiVertexInfo & left, const Simplex_handle & right) ->bool
					{
						return left.hDelaunaySimplex<right;
					});
					if (NewFaceIterator == vInfo.end())
						continue;
					if (NewFaceIterator->hDelaunaySimplex == hNeighbor)
					{
						//debug temp
						//std::cerr<<iter->VerticesVectorLocation<<" \t"<<(d+1)<<" \t"<<NeighborCount<<'\n';
						pLinkage->at(iter->VerticesVectorLocation*(d + 1) + NeighborCount) = NewFaceIterator->VerticesVectorLocation;
						NeighborCount++;
					}
				}
				//for(; NeighborCount<d+1; NeighborCount++)
				//	pLinkage->at( iter->VerticesVectorLocation*(d+1) + NeighborCount ) = -1;
			}
		}
	}
	//debug temp
	//std::cout << "Num iteration = " << NumIteration << '\n';
}

const double Perturbation = 0;
double PolygonVolume_ConvexHull(const std::vector<GeometryVector> & Vertices, DimensionType d, RandomGenerator & gen)
{
	if (Vertices.size() == 0)
		return 0;
	ConvexHull De(d);
	GeometryVector center = SameDimensionZeroVector(Vertices[0]);
	for (auto iter = Vertices.begin(); iter != Vertices.end(); iter++)
	{
		GeometryVector temp(*iter);
		for (DimensionType i = 0; i<d; i++)
			temp.x[i] += Perturbation * temp.x[i] * (gen.RandomDouble() - 0.5);
		center.AddFrom(temp);
		cPoint p(d, temp.x, temp.x + d);
		De.insert(p);
	}
	center = center * (1.0 / Vertices.size());
	double result = 0;
	std::list<ConvexHull::Facet_handle> facets = De.all_facets();
	for (auto iter = facets.begin(); iter != facets.end(); iter++)
	{
		std::vector<GeometryVector> vecs;
		for (DimensionType i = 0; i<d; i++)
		{
			GeometryVector vi = CGALPoint2GeometryVector((De.point_of_facet(*iter, i)), d);
			vecs.push_back(vi - center);
		}
		result += ::SimplexVolume(&vecs[0], d);
	}
	return result;
}
double PolygonVolume(const std::vector<GeometryVector> & Vertices, DimensionType d)
{
	if (d == 1)
	{
		std::vector<double> coords;
		for (auto iter = Vertices.begin(); iter != Vertices.end(); iter++)
			coords.push_back(iter->x[0]);
		//std::nth_element(coords.begin(), coords.begin(), coords.end());
		//std::nth_element(coords.begin()+1, coords.end()-1, coords.end());
		//return coords.back()-coords.front();
		return *std::max_element(coords.begin(), coords.end()) - *std::min_element(coords.begin(), coords.end());
	}
	else
	{
		RandomGenerator gen(1);
		return PolygonVolume_ConvexHull(Vertices, d, gen);
	}
}



void PlotVoronoi2D(std::string prefix, const Configuration & List, const char * edgeColor, const std::string & remark)
{
	const Configuration * list = &List;
#ifdef USE_MATHGL

	if (list->GetDimension() == 2)
	{
		//calculate the range of the plot
		GeometryVector a0(0, 0);
		GeometryVector a1 = a0 + list->GetBasisVector(0);
		GeometryVector a2 = a0 + list->GetBasisVector(1);
		GeometryVector a3 = a2 + list->GetBasisVector(0);
		double min = 0, max = 0;
		if (min>a0.x[0])
			min = a0.x[0];
		if (min>a0.x[1])
			min = a0.x[1];
		if (max<a0.x[0])
			max = a0.x[0];
		if (max<a0.x[1])
			max = a0.x[1];
		if (min>a1.x[0])
			min = a1.x[0];
		if (min>a1.x[1])
			min = a1.x[1];
		if (max<a1.x[0])
			max = a1.x[0];
		if (max<a1.x[1])
			max = a1.x[1];
		if (min>a2.x[0])
			min = a2.x[0];
		if (min>a2.x[1])
			min = a2.x[1];
		if (max<a2.x[0])
			max = a2.x[0];
		if (max<a2.x[1])
			max = a2.x[1];
		if (min>a3.x[0])
			min = a3.x[0];
		if (min>a3.x[1])
			min = a3.x[1];
		if (max<a3.x[0])
			max = a3.x[0];
		if (max<a3.x[1])
			max = a3.x[1];

		//prepare the graph
		mglGraph gr;
		double delta = max - min;
		min = min - 0.1*delta;
		max = max + 0.1*delta;
		gr.SetRanges(min, max, min, max);
		gr.Aspect(1.0, 1.0);
		gr.SetSize(800, 800);
		//gr.SetOrigin(0, 0);
		gr.Puts(mglPoint(0.5*(max + min), min + 1.1*(max - min)), remark.c_str());

		//plot crystal BasisVectors
		{
			mglData x(2), y(2);
			x.a[0] = a0.x[0];
			y.a[0] = a0.x[1];
			x.a[1] = a1.x[0];
			y.a[1] = a1.x[1];
			gr.Plot(x, y, "k-");
			x.a[0] = a1.x[0];
			y.a[0] = a1.x[1];
			x.a[1] = a3.x[0];
			y.a[1] = a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0] = a2.x[0];
			y.a[0] = a2.x[1];
			x.a[1] = a3.x[0];
			y.a[1] = a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0] = a2.x[0];
			y.a[0] = a2.x[1];
			x.a[1] = a0.x[0];
			y.a[1] = a0.x[1];
			gr.Plot(x, y, "k-");
		}
		//plot particles
		mglData x(list->NumParticle()), y(list->NumParticle());
		std::vector<double> xx, yy;
		for (int i = 0; i<list->NumParticle(); i++)
		{
			GeometryVector now = list->GetCartesianCoordinates(i);
			x.a[i] = now.x[0];
			y.a[i] = now.x[1];

			std::vector<GeometryVector> vertex;
			std::vector<signed long> linkage;
			::GetVoronoiCell(List, i, vertex, &linkage);
			for (size_t i = 0; i<linkage.size(); i++)
			{
				if (linkage[i] == -1)
					continue;
				size_t Index1 = linkage[i];
				size_t Index2 = i / 3;
				if (Index1<Index2)
					continue;
				mglData x(2), y(2);
				x.a[0] = vertex[Index1].x[0];
				y.a[0] = vertex[Index1].x[1];
				x.a[1] = vertex[Index2].x[0];
				y.a[1] = vertex[Index2].x[1];
				gr.Plot(x, y, edgeColor);
			}
			for (auto iter = vertex.begin(); iter != vertex.end(); iter++)
			{
				xx.push_back(iter->x[0]);
				yy.push_back(iter->x[1]);
			}
		}
		gr.Plot(x, y, "b .");
		mglData xxx(&xx[0], xx.size()), yyy(&yy[0], yy.size());
		//gr.Plot(xxx, yyy, "r .");

		prefix += std::string(FigureFormat);
		gr.WriteFrame(prefix.c_str());
		if (AlsoWriteEPS)
		{
			prefix.resize(prefix.size() - 4);
			prefix += std::string(".eps");
			gr.WriteEPS(prefix.c_str());
		}
	}
	else
	{
		std::cerr << "error in PlotVoronoi2D : unsupported dimension\n";
}
#endif
}

//PlotVoronoi2D+PlotPacking
void PlotVoronoi2D(std::string prefix, const SpherePacking & List, const char * sphereColor, const char * edgeColor, const std::string & remark)
{
	const SpherePacking * list = &List;
#ifdef USE_MATHGL

	if (list->GetDimension() == 2)
	{
		//calculate the range of the plot
		GeometryVector a0(0, 0);
		GeometryVector a1 = a0 + list->GetBasisVector(0);
		GeometryVector a2 = a0 + list->GetBasisVector(1);
		GeometryVector a3 = a2 + list->GetBasisVector(0);
		double min = 0, max = 0;
		if (min>a0.x[0])
			min = a0.x[0];
		if (min>a0.x[1])
			min = a0.x[1];
		if (max<a0.x[0])
			max = a0.x[0];
		if (max<a0.x[1])
			max = a0.x[1];
		if (min>a1.x[0])
			min = a1.x[0];
		if (min>a1.x[1])
			min = a1.x[1];
		if (max<a1.x[0])
			max = a1.x[0];
		if (max<a1.x[1])
			max = a1.x[1];
		if (min>a2.x[0])
			min = a2.x[0];
		if (min>a2.x[1])
			min = a2.x[1];
		if (max<a2.x[0])
			max = a2.x[0];
		if (max<a2.x[1])
			max = a2.x[1];
		if (min>a3.x[0])
			min = a3.x[0];
		if (min>a3.x[1])
			min = a3.x[1];
		if (max<a3.x[0])
			max = a3.x[0];
		if (max<a3.x[1])
			max = a3.x[1];

		//prepare the graph
		mglGraph gr;
		double delta = max - min;
		min = min - 0.1*delta;
		max = max + 0.1*delta;
		gr.SetRanges(min, max, min, max);
		gr.Aspect(1.0, 1.0);
		gr.SetSize(800, 800);
		//gr.SetOrigin(0, 0);
		gr.Puts(mglPoint(0.5*(max + min), min + 1.1*(max - min)), remark.c_str());

		//plot crystal BasisVectors
		{
			mglData x(2), y(2);
			x.a[0] = a0.x[0];
			y.a[0] = a0.x[1];
			x.a[1] = a1.x[0];
			y.a[1] = a1.x[1];
			gr.Plot(x, y, "k-");
			x.a[0] = a1.x[0];
			y.a[0] = a1.x[1];
			x.a[1] = a3.x[0];
			y.a[1] = a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0] = a2.x[0];
			y.a[0] = a2.x[1];
			x.a[1] = a3.x[0];
			y.a[1] = a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0] = a2.x[0];
			y.a[0] = a2.x[1];
			x.a[1] = a0.x[0];
			y.a[1] = a0.x[1];
			gr.Plot(x, y, "k-");
		}
		//plot particles
		mglData x(list->NumParticle()), y(list->NumParticle());
		std::vector<double> xx, yy;
		for (int i = 0; i < list->NumParticle(); i++)
		{

			double r = list->GetCharacteristics(i); //^M
				GeometryVector now = list->GetCartesianCoordinates(i); //^M
				::PlotDisk2D(gr, now.x[0], now.x[1], r, sphereColor);

			x.a[i] = now.x[0];
			y.a[i] = now.x[1];

			std::vector<GeometryVector> vertex;
			std::vector<signed long> linkage;
			::GetVoronoiCell(List, i, vertex, &linkage);
			//DRAW CELL EDGES
			for (size_t i = 0; i<linkage.size(); i++)
			{
				if (linkage[i] == -1)
					continue;
				size_t Index1 = linkage[i];
				size_t Index2 = i / 3;
				if (Index1<Index2)
					continue;
				mglData x(2), y(2);
				x.a[0] = vertex[Index1].x[0];
				y.a[0] = vertex[Index1].x[1];
				x.a[1] = vertex[Index2].x[0];
				y.a[1] = vertex[Index2].x[1];
				gr.Plot(x, y, edgeColor);
			}
			for (auto iter = vertex.begin(); iter != vertex.end(); iter++)
			{
				xx.push_back(iter->x[0]);
				yy.push_back(iter->x[1]);
			}
		}
		gr.Plot(x, y, "k .");
		{//draw vertices
		 //mglData xxx(&xx[0], xx.size()), yyy(&yy[0], yy.size());
		 //gr.Plot(xxx, yyy, "r .");
		}
		prefix += std::string(FigureFormat);
		gr.WriteBMP(prefix.c_str());
		if (AlsoWriteEPS)
		{
			prefix.resize(prefix.size() - 4);
			prefix += std::string(".eps");
			gr.WriteEPS(prefix.c_str());
		}
	}
	else
	{
		std::cerr << "error in PlotVoronoi2D : unsupported dimension\n";
	}
#endif
}
		


void DisplayConfiguration_3D_WithVoronoi(const Configuration & list, double Radius, long Num)
{
	assert(list.GetDimension() == 3);
	long & TargetI = Num;
	//debug temp
	//std::cin>>TargetI;
	double typicalLength = std::pow(list.PeriodicVolume() / list.NumParticle(), 1.0 / 3.0);
	for (size_t i = 0; i<list.NumParticle(); i++)
	{
		//Configuration::particle * pa=list.GetParticle(i);
		list.IterateThroughNeighbors(list.GetRelativeCoordinates(i), typicalLength, [&typicalLength, &i](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
		{
			if (SourceAtom != i)
			{
				double l = std::sqrt(shift.Modulus2());
				if (l<typicalLength)
					typicalLength = l;
			}
		});
	}

	std::vector<sphere> spheres;
	std::vector<line> lines;

	std::vector<GeometryVector> a(8, GeometryVector(0.0, 0.0, 0.0));
	a[1] = a[0] + list.GetBasisVector(0);
	a[2] = a[0] + list.GetBasisVector(1);
	a[3] = a[0] + list.GetBasisVector(2);
	a[4] = a[1] + list.GetBasisVector(1);
	a[5] = a[1] + list.GetBasisVector(2);
	a[6] = a[2] + list.GetBasisVector(2);
	a[7] = a[4] + list.GetBasisVector(2);
	GeometryVector shift = a[7] * 0.5;
	double one_over_scale = 1.0 / std::sqrt(shift.Modulus2());
	for (int i = 0; i<8; i++)
	{
		a[i].MinusFrom(shift);
		a[i].MultiplyFrom(one_over_scale);
	}
	int point1[12] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 5, 6 };
	int point2[12] = { 1, 2, 3, 4, 5, 4, 6, 5, 6, 7, 7, 7 };
	for (int i = 0; i<12; i++)
	{
		line temp;
		temp.x1 = a[point1[i]].x[0];
		temp.y1 = a[point1[i]].x[1];
		temp.z1 = a[point1[i]].x[2];
		temp.x2 = a[point2[i]].x[0];
		temp.y2 = a[point2[i]].x[1];
		temp.z2 = a[point2[i]].x[2];
		temp.blue = 0.0;
		temp.red = 0.0;
		temp.green = 0.0;
		temp.transparency = 1.0;
		temp.width = 0.06*typicalLength*one_over_scale;
		lines.push_back(temp);
	}

	for (size_t i = 0; i<list.NumParticle(); i++)
	{
		//auto pA=list.GetParticle(i);
		GeometryVector loc = (list.GetCartesianCoordinates(i) - shift)*one_over_scale;
		sphere temp;
		temp.x = loc.x[0];
		temp.y = loc.x[1];
		temp.z = loc.x[2];

		char chr1 = list.GetCharacteristics(i).name[0];
		if (chr1 >= 'A' && chr1 <= 'Z')
		{
			temp.red = ((chr1 - 'A') / 6)*0.1 + 0.3;
			temp.green = ((chr1 - 'A') % 6)*0.1 + 0.3;
			char chr2 = list.GetCharacteristics(i).name[1];
			if (chr2 >= 'a' && chr2 <= 'z')
				temp.blue = (chr2 - 'a')*0.02 + 0.4;
			else
				temp.blue = 0.35;
		}
		else
		{
			temp.red = 0.25;
			temp.green = 0.25;
			temp.blue = 0.25;
		}
		temp.transparency = 0.97;

		if (Radius == 0.0)
			temp.radius = 0.3*typicalLength*one_over_scale;
		else
			temp.radius = Radius * one_over_scale;

		spheres.push_back(temp);

		if (i == TargetI || TargetI == -1)
		{
			std::vector<GeometryVector> vertex;
			std::vector<signed long> linkage;
			::GetVoronoiCell(list, i, vertex, &linkage);
			for (size_t i = 0; i<linkage.size(); i++)
			{
				if (linkage[i] == -1)
					continue;
				size_t Index1 = linkage[i];
				size_t Index2 = i / 4;
				if (Index1<Index2)
					continue;
				GeometryVector end1 = (vertex[Index1] - shift)*one_over_scale;
				GeometryVector end2 = (vertex[Index2] - shift)*one_over_scale;
				line temp;
				temp.x1 = end1.x[0];
				temp.x2 = end2.x[0];
				temp.y1 = end1.x[1];
				temp.y2 = end2.x[1];
				temp.z1 = end1.x[2];
				temp.z2 = end2.x[2];
				temp.blue = 1.0;
				temp.red = 0.0;
				temp.green = 0.0;
				temp.transparency = 1.0;
				temp.width = 0.03*typicalLength*one_over_scale;
				lines.push_back(temp);
			}
			for (size_t i = 0; i<vertex.size(); i++)
			{
				GeometryVector loc = (vertex[i] - shift)*one_over_scale;
				sphere temp;
				temp.x = loc.x[0];
				temp.y = loc.x[1];
				temp.z = loc.x[2];

				temp.blue = 1.0;
				temp.red = 0.0;
				temp.green = 0.0;
				temp.transparency = 0.97;

				if (Radius == 0.0)
					temp.radius = 0.1*typicalLength*one_over_scale;
				else
					temp.radius = Radius * one_over_scale*0.33;

				spheres.push_back(temp);
			}
		}

	}

	DisplaySpheres(spheres, lines);
}
#include <omp.h>
void VoronoiVolumeDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, size_t SampleSize, double ResolutionPreference)
{
	//get the pair distances list
	if (Verbosity>2)
		std::cout << "Generate VoronoiVolumeDistrubution";
	progress_display pd(NumConfigs);

	std::vector<double> distances;

	DimensionType dim;

	size_t NumConfigsSampled = 0;

	omp_lock_t lock;
	omp_init_lock(&lock);

	auto GetVoronoiSizeFunction = [](const Configuration & c, size_t i, DimensionType dim) ->double
	{
		std::vector<GeometryVector> vertices;
		::GetVoronoiCell(c, i, vertices);
		return PolygonVolume(vertices, dim);
	};
	for (size_t j = 0; j<NumConfigs; j++)
	{
		pd++;
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration Config = GetConfigsFunction(j);
		if (j == 0)
			dim = Config.GetDimension();
		size_t NumParticle = Config.NumParticle();
		double TypicalLength = std::pow(Config.PeriodicVolume() / Config.NumParticle(), 1.0 / Config.GetDimension())*1.5;

		long end = Config.NumParticle();
#pragma omp parallel 
		{
			Configuration ThreadConfig(Config);
#pragma omp for schedule(dynamic)
			for (long j = 0; j<end; j++)
			{
				if (Verbosity>4)
					std::cout << "Doing Particle" << j << '\n';
				double result = GetVoronoiSizeFunction(ThreadConfig, j, dim);
				omp_set_lock(&lock);
				distances.push_back(result);
				omp_unset_lock(&lock);
			}
		}
		NumConfigsSampled++;

		if (distances.size()>SampleSize)
			break;
	}
	omp_destroy_lock(&lock);

	//if(Verbosity>2)
	//	std::cout<<"generating bins...";
	//prepare for the iteration
	Result.clear();
	//iteration
	HistogramGenerator HGen;
	HGen.GenerateBins(distances, static_cast<double>(NumConfigs) / NumConfigsSampled, ResolutionPreference);

	//Re-do the bin counts
	//if(Verbosity>2)
	//	std::cout<<"counting pair distances...";

	for (size_t j = NumConfigsSampled; j<NumConfigs; j++)
	{
		pd++;
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration Config = GetConfigsFunction(j);
		size_t NumParticle = Config.NumParticle();

		long end = Config.NumParticle();
#pragma omp parallel 
		{
			Configuration ThreadConfig(Config);
#pragma omp for schedule(dynamic)
			for (long j = 0; j<end; j++)
			{
				if (Verbosity>4)
					std::cout << "Doing Particle" << j << '\n';
				HGen.Report(GetVoronoiSizeFunction(ThreadConfig, j, dim));
			}
		}
	}

	size_t TotalCount = 0;
	for (auto iter = HGen.bins.begin(); iter != HGen.bins.end(); iter++)
		TotalCount += iter->Count;

	for (auto iter = HGen.bins.begin(); iter != HGen.bins.end(); iter++)
	{
		GeometryVector temp(4);
		if (iter != HGen.bins.end() - 1)
			temp.x[0] = ((iter[0].Start) + (iter[1].Start)) / 2.0;//r 
		else
			continue;

		if (iter == HGen.bins.begin())
		{
			if (iter->Start == 0.0 && (iter[2].Start - iter[1].Start)>iter[1].Start && iter->Count == 0)
				continue;//there shouldn't be a bin at here
		}
		temp.x[2] = (temp.x[0] - iter[0].Start);//\delta r

		temp.x[1] = static_cast<double>(iter->Count) / TotalCount / (2 * temp.x[2]);//p(r)

		if (iter->Count != 0)
			temp.x[3] = temp.x[1] / std::sqrt(static_cast<double>(iter->Count));//\delta p(r)
		else
			temp.x[3] = 0.0;
		Result.push_back(temp);
	}

	if (Verbosity>2)
		std::cout << "done!\n";
}
void VoronoiNumSidesDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig, std::vector<GeometryVector> & Result)
{
	//get the pair distances list
	if (Verbosity>2)
		std::cout << "Generate VoronoiNumSidesDistrubution";
	std::map<size_t, size_t> sides;
	progress_display pd(NumConfig);
	omp_lock_t pdlock;
	omp_init_lock(&pdlock);
	omp_lock_t sideslock;
	omp_init_lock(&sideslock);
	size_t SumNumParticle = 0;
#pragma omp parallel for
	for (signed long i = 0; i < NumConfig; i++)
	{
		Configuration c = GetConfigsFunction(i);
		size_t np = c.NumParticle();
#pragma omp atomic
		SumNumParticle += np;
		for (size_t j = 0; j < np; j++)
		{
			std::vector<GeometryVector> vertices;
			std::vector<signed long> linkage;
			GetVoronoiCell(c, j, vertices, &linkage);
			size_t nside = 0;
			for (auto iter = linkage.begin(); iter != linkage.end(); ++iter)
				if (*iter != -1)
					nside++;
			assert(nside % 2 == 0);
			nside /= 2;
			omp_set_lock(&sideslock);
			std::map<size_t, size_t>::iterator iter = sides.find(nside);
			if (iter == sides.end())
				sides.insert(std::make_pair(nside, (size_t)(1)));
			else
				iter->second++;
			omp_unset_lock(&sideslock);
		}
		omp_set_lock(&pdlock);
		pd++;
		omp_unset_lock(&pdlock);
	}
	omp_destroy_lock(&sideslock);
	omp_destroy_lock(&pdlock);
	for (auto iter = sides.begin(); iter != sides.end(); iter++)
		Result.push_back(GeometryVector(iter->first, (double)(iter->second) / SumNumParticle));

	if (Verbosity>2)
		std::cout << "done!\n";
}


//C(r)=<(V(x)-<V>)(V(x+r)-<V>)>_x
void VoronoiVolumeCorrelation(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, double MaxDistance, double RPrecision)
{
	//get the pair distances list
	if (Verbosity>2)
		std::cout << "Generate VoronoiVolumeCorrelation";
	progress_display pd(NumConfigs);

	std::vector<double> distances;

	DimensionType dim;

	auto GetVoronoiSizeFunction = [](const Configuration & c, size_t i, DimensionType dim) ->double
	{
		std::vector<GeometryVector> vertices;
		::GetVoronoiCell(c, i, vertices);
		return PolygonVolume(vertices, dim);
	};
	//double SumVolume=0;
	//size_t SumNumParticle=0;

	//this function calculates <(V(x)-<V>)(V(x+r)-<V>)>_x
	//this struct is a single data of (V(x)-<V>)(V(x+r)-<V>)
	struct DataType
	{
		double r;
		double vv;
	};
	std::vector<DataType> vData;

	for (size_t j = 0; j<NumConfigs; j++)
	{
		pd++;
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration Config = GetConfigsFunction(j);
		if (j == 0)
			dim = Config.GetDimension();
		size_t NumParticle = Config.NumParticle();
		//SumNumParticle+=NumParticle;
		double AverageVolume = Config.PeriodicVolume() / Config.NumParticle();
		std::vector<double> VVolume(NumParticle, 0.0);

		long end = NumParticle;
#pragma omp parallel 
		{
			Configuration ThreadConfig(Config);
#pragma omp for schedule(dynamic)
			for (long j = 0; j<end; j++)
			{
				if (Verbosity>4)
					std::cout << "Doing Particle" << j << '\n';
				double result = GetVoronoiSizeFunction(ThreadConfig, j, dim);
				VVolume[j] = result - AverageVolume;
				//#pragma omp atomic
				//				SumVolume+=result;
			}
		}

		for (long j = 0; j<end; j++)
		{
			auto func = [&](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
			{
				DataType temp;
				if (SourceAtom == j && shift.Modulus2()<1e-10)
					temp.r = 0.0;
				else
					temp.r = std::sqrt(shift.Modulus2());
				temp.vv = VVolume[j] * VVolume[SourceAtom];
				vData.push_back(temp);
			};
			GeometryVector Rel = Config.GetRelativeCoordinates(j);
			Config.IterateThroughNeighbors(Rel, MaxDistance, func);
		}
	}

	//	double AverageVolume=SumVolume/SumNumParticle;
	auto CompareFunc = [](const DataType & left, const DataType & right) ->bool
	{
		return left.r<right.r;
	};
	std::sort(vData.begin(), vData.end(), CompareFunc);

	size_t SameRCount = 0;
	double Sum = 0.0;
	double Sum2 = 0.0;
	double CurrentR = 0.0;
	double SumR2 = 0.0;
	//if(KPrecision==0.0)
	//	KPrecision = LengthPrecision/tempList.GetBasisVector(0).Modulus2();
	for (auto iter = vData.begin(); iter != vData.end(); iter++)
	{
		double c = iter->vv;
		double c2 = c * c;
		double r = iter->r;
		if ((r - CurrentR < RPrecision&& CurrentR>0.0) || (r == 0.0 && CurrentR == 0.0))
		{
			SameRCount++;
			Sum += c;
			Sum2 += c2;
			SumR2 += r * r;
		}
		else
		{
			if (SameRCount>0)
			{
				GeometryVector temp(4);
				temp.x[0] = std::sqrt(SumR2 / SameRCount);
				temp.x[1] = Sum / SameRCount;
				if (CurrentR>0.0)
					temp.x[2] = RPrecision;
				else
					temp.x[2] = 0;
				temp.x[3] = std::sqrt(Sum2 / (SameRCount)-temp.x[1] * temp.x[1]) / (SameRCount);
				Result.push_back(temp);
			}

			Sum = c;
			Sum2 = c * c;
			SameRCount = 1;
			CurrentR = r;
			SumR2 = r * r;
		}
	}
	if (SameRCount>0)
	{
		GeometryVector temp(4);
		temp.x[0] = std::sqrt(SumR2 / SameRCount);
		temp.x[1] = Sum / SameRCount;
		temp.x[2] = RPrecision;
		temp.x[3] = std::sqrt(Sum2 / (SameRCount)-temp.x[1] * temp.x[1]) / (SameRCount);
		Result.push_back(temp);
	}

	if (Verbosity>2)
		std::cout << "done!\n";
}

#else

void GetVoronoiCell(const PeriodicCellList<Empty> & c, size_t i, std::vector<GeometryVector> & vertices, std::vector<signed long> * pLinkage, std::set<size_t> * pNeighbors)
{
	std::cerr << "Error : CGAL not enabled!\n";
	return;
}
double PolygonVolume(std::vector<GeometryVector> Vertices, DimensionType d)
{
	std::cerr << "Error : CGAL not enabled!\n";
	return 0;
}
void PlotVoronoi2D(std::string prefix, const Configuration & List, const std::string & remark)
{
	std::cerr << "Error : CGAL not enabled!\n";
	return;
}

//Num=-1 : display voronoi of all particles Num>=0 : display voronoi of a specific particle
void DisplayConfiguration_3D_WithVoronoi(const Configuration & list, double Radius, long Num)
{
	std::cerr << "Error : CGAL not enabled!\n";
	return;
}

void VoronoiVolumeDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, size_t SampleSize, double ResolutionPreference)
{
	std::cerr << "Error : CGAL not enabled!\n";
	return;
}
void VoronoiVolumeCorrelation(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, double MaxDistance, double RPrecision)
{
	std::cerr << "Error : CGAL not enabled!\n";
	return;
}

void VoronoiNumSidesDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig, std::vector<GeometryVector> & Result)
{
	std::cerr << "Error : CGAL not enabled!\n";
	return;
}
#endif


//param[in] input	initial configuration
//param[in] output	an empty configuration with rotated basis vectors in the following way 
//					b[0].x[1] = b[0].x[2] = 0	b[1].x[2] = 0;
//param[out]		true: basis vectors are changed!	false: basis vectors aren't changed.
bool transformBasis(const PeriodicCellList<Empty> & input, Configuration & output) {
	double tolerance = 1e-20; //numerical tolerance!!



	DimensionType d = input.GetDimension();
	//Initial basis vectors!!
	std::vector<GeometryVector> b (d);
	for (int i = 0; i < (int)d; i++) {
		b[i] = input.GetBasisVector(i);
	}

	bool change;
	if (d == 1) {
		return false;
	}
	else if (d == 2) {

		change = (fabs(b[0].x[1]) > tolerance);
		if (change) {//Rotate basis vectors!!
			double bx = sqrt(b[0].Modulus2());
			double bxy = b[0].Dot(b[1]) / sqrt(b[0].Modulus2()), byy = sqrt(b[1].Modulus2() - bxy * bxy);

			//Final basis vectors!!
			std::vector<GeometryVector> b_(d);

			b_[0].Dimension = d;
			b_[0].x[0] = bx; b_[0].x[1] = 0.0;
			b_[1].Dimension = d;
			b_[1].x[0] = bxy; b_[1].x[1] = byy;
			output = Configuration(d, &b_[0], 1.);

			{//Debug
				for (int j = 0; j< d; j++)
					std::cout << "|b_[" << j << "]| = " << sqrt(b_[j].Modulus2()) << "\n";
			}
		}
		return change;
	}
	else if (d == 3) {

		change = (fabs(b[0].x[1]) > tolerance) || (fabs(b[0].x[2]) > tolerance) || (fabs(b[1].x[2]) > tolerance);
		if (change) {
			double bx = sqrt(b[0].Modulus2());
			double bxy = b[0].Dot(b[1]) / sqrt(b[0].Modulus2()), byy = sqrt(b[1].Modulus2() - bxy * bxy);

			double z_x = 0.0, z_y = 0.0, z_z = 1.0; //strip...
													//Final basis vectors!!
			std::vector<GeometryVector> b_(d, GeometryVector(0.0, 0.0, 0.0));

			b_[0].x[0] = bx;
			b_[1].x[0] = bxy; b_[1].x[1] = byy;
			b_[2] = b[0].Dot(b[2]) / b[0].Modulus2() * b_[0] + b[1].Dot(b[2]) / b[1].Modulus2() *b_[1];

			GeometryVector cross(b[0].x[1] * b[1].x[2] - b[0].x[2] * b[1].x[1], b[0].x[2] * b[1].x[0] - b[0].x[0] * b[1].x[2], b[0].x[0] * b[1].x[1] - b[0].x[1] * b[1].x[0]);
			b_[2].x[2] = cross.Dot(b[2]) / sqrt(cross.Modulus2());
			{//Debug
				for (int j = 0; j < d; j++) {
					std::cout << "|b_[" << j << "]| = " << sqrt(b_[j].Modulus2()) << "\n";
					std::cout << b_[j];
				}
			}
			output = Configuration(d, &b_[0], 1.);
		}

		return change;

	}
	else {
		std::cout << "transformBasis::Dimension is too high! \n";
		return false;
	}
}




#if defined(USE_VORO)
void container_leak::GetVORONOICELL(int i, std::vector<GeometryVector> & v, std::vector<signed long> * pLinkage, std::set<size_t> * pNeighbors) {
	v.resize(0);
	
	int ijk = this->id_leak[i].blockId;
	int q = this->id_leak[i].order;
	int d = this->c_.GetDimension();

	std::vector<double> vList;	//A list of vertex coordinates
	GeometryVector vi(static_cast<int>(d));
	GeometryVector x0(vi);	//Particle position
	
	/*double *pp = this->p[ijk] + this->ps*q;
	if (d == 2) {
		x0.x[0] = *pp;	x0.x[1] = pp[1];
	}
	else if (d == 3) {
		x0.x[0] = *pp;	x0.x[1] = pp[1];	x0.x[2] = pp[2];
	}*/
	if (d == 2) {
		x0.x[0] = this->p[ijk][this->ps*q];	x0.x[1] = this->p[ijk][this->ps*q+1];
	}
	else if (d == 3) {
		x0.x[0] = this->p[ijk][this->ps*q];	x0.x[1] = this->p[ijk][this->ps*q+1];	x0.x[2] = this->p[ijk][this->ps*q+2];
	}

	if (this->transform == true)//Do inverse coordinate transform!
		x0 = this->c0.RelativeCoord2CartesianCoord(c_.CartesianCoord2RelativeCoord(x0));
	
	if (pLinkage == nullptr) {
		voro::voronoicell vcell;
		if (this->compute_cell(vcell, ijk, q)) {
			vcell.vertices(vList);
			for (int j = 0; j < vcell.p; j++) {
				if (d == 2 && vList[3 * j + 2] > 0.0) {
					vi.x[0] = vList[3 * j];	vi.x[1] = vList[3 * j + 1];				
				}
				else if (d == 3) {
					vi.x[0] = vList[3 * j];	vi.x[1] = vList[3 * j + 1];	vi.x[2] = vList[3 * j + 2];
				}

				if ((d == 2 && vList[3 * j + 2]>0.0) || d==3) {
					if (this->transform)//Do inverse coordinate transform!
						vi = this->c0.RelativeCoord2CartesianCoord(c_.CartesianCoord2RelativeCoord(vi));

					vi = vi + x0;
					v.push_back(vi);
				}
			}
		}
		else {
			std::cout << "Something wrong!\n";
			return;
		}

	}
}

container_leak::container_leak(const PeriodicCellList<Empty> & cInit, const std::vector<GeometryVector> & newBasis, bool trans, int nx, int ny,int nz, int mem):
	container_leak::container_leak(newBasis[0].x[0], newBasis[1].x[0], newBasis[1].x[1], newBasis[2].x[0], newBasis[2].x[1], newBasis[2].x[2], nx, ny, nz, mem)
{
	//Debug
	//int count = 0;
	int d = cInit.GetDimension();
	this->c0 = Configuration(cInit, "a");	c0.RemoveParticles();
	this->transform = trans;
	this->id_leak = std::vector<IDs>(cInit.NumParticle(), IDs());
	{
		std::vector<GeometryVector> basis;
			if (d == 2) {
				basis.emplace_back(GeometryVector(newBasis[0].x[0], newBasis[0].x[1]));
				basis.emplace_back(GeometryVector(newBasis[1].x[0], newBasis[1].x[1]));
			}
			else if (d == 3) {
				basis.emplace_back(GeometryVector(newBasis[0].x[0], newBasis[0].x[1], newBasis[0].x[2]));
				basis.emplace_back(GeometryVector(newBasis[1].x[0], newBasis[1].x[1], newBasis[1].x[2]));
				basis.emplace_back(GeometryVector(newBasis[2].x[0], newBasis[2].x[1], newBasis[2].x[2]));
			}
		this->c_ = Configuration(d, &basis[0], sqrt(basis[0].Modulus2()));
	}
	if (d == 2 || d == 3) {
		//Input particles
		int ijk; double x, y, z;	GeometryVector temp;
		for (int i = 0; i < cInit.NumParticle(); i++) {
			if (trans) { temp = c_.RelativeCoord2CartesianCoord(cInit.GetRelativeCoordinates(i)); }
			else { temp = cInit.GetCartesianCoordinates(i); }

			double x = temp.x[0], y = temp.x[1], z = 0.5;
			if (d == 3)
				z = temp.x[2];

			ijk = this->put_ijk(i, x, y, z);
			this->id_leak[i].blockId = ijk;
			this->id_leak[i].order = (this->co[ijk]);
			//Debug
			//if(this->id_leak[i].order >= mem){
			//	count++;
			//}
		}
		for (int i = 0; i < id_leak.size(); i++)
			this->id_leak[i].order--;

		//std::cout << "overflow = ?" << count << "\n";
	}
}

container_leak::container_leak(const container_leak & source):
	container_leak::container_leak(source.bx, source.bxy, source.by, source.bxz, source.byz, source.bz, source.nx, source.ny, source.nz, source.init_mem) 
{
	this->c0 = Configuration(source.c0);
	this->c_ = Configuration(source.c_);
	this->transform = source.transform;
	this->id_leak = std::vector<IDs>(source.id_leak);

	double x, y, z;	int n;
	voro::c_loop_all_periodic routine(source);
	
	if (routine.start() && source.id_leak.size() > 0) {
		do {
			n = routine.pid();
			x = routine.x();
			y = routine.y();
			z = routine.z();
			this->put(n, x, y, z);
		} while (routine.inc());
	}
}

void GetConversionParam(const PeriodicCellList<Empty> &c, int approxNumPrts, std::vector<GeometryVector> & basis, bool &trans, int & nx, int & ny, int &nz, int mem) {
	Configuration c_new;
	int d = c.GetDimension();
	trans = transformBasis(c, c_new);
	basis = std::vector<GeometryVector>(3, GeometryVector(0.0, 0.0, 0.0));
	for (int i = 0; i < 3; i++) {
		basis[i].x[i] = 1.0;
		if (i < d) {
			for (int j = 0; j < d; j++) {
				if (trans)
					basis[i].x[j] = c_new.GetBasisVector(i).x[j];
				else
					basis[i].x[j] = c.GetBasisVector(i).x[j];
			}
		}
	}

	nx = (int)ceil(pow(approxNumPrts / mem, 1.0 / d));
	ny = nx;
	if (d == 2)
		nz = 1;
	else if (d == 3)
		nz = nx;
}

void GetConversionParam(const PeriodicCellList<Empty> &c, std::vector<GeometryVector> & basis, bool &trans, int & nx, int & ny, int &nz, int mem) {
	GetConversionParam(c, c.NumParticle(), basis, trans, nx, ny, nz, mem);
}




container_leak Config2Container(const PeriodicCellList<Empty> & c) {
	int d = c.GetDimension();
	int mem = 10;
	Configuration c_, c0;
	c0 = Configuration(c, "a");	c0.RemoveParticles();

	bool trans = transformBasis(c, c_);
	if (trans == false) {
		c_ = Configuration(c0);
	}

	if (d == 2 || d == 3) {
		int n = (int)ceil(pow(c.NumParticle() / (double)mem, 1.0 / (double)d) );
		int nz = (d == 2) ? 1 : n;
		double z1, z2, z3;
		if (d == 2) {
			z1 = 0;	z2 = 0;	z3 = 1.0;
		}
		else {
			z1 = c_.GetBasisVector(2).x[0];	z2 = c_.GetBasisVector(2).x[1];	z3 = c_.GetBasisVector(2).x[2];
		}
		container_leak unitCell(c_.GetBasisVector(0).x[0], c_.GetBasisVector(1).x[0], c_.GetBasisVector(1).x[1], z1, z2, z3, n, n, nz, mem);

		unitCell.c_ = Configuration(c_);
		unitCell.transform = trans;
		unitCell.id_leak = std::vector<IDs>(c.NumParticle(), IDs());
		//Input particles
		int ijk; double x, y, z;	GeometryVector temp;
		for (int i = 0; i < c.NumParticle(); i++) {
			if (trans) { temp = c_.RelativeCoord2CartesianCoord(c.GetRelativeCoordinates(i)); }
			else { temp = c.GetCartesianCoordinates(i); }

			double x = temp.x[0], y = temp.x[1], z = 0.5;
			if (d == 3)
				z = temp.x[2];

			ijk = unitCell.put_ijk(i, x, y, z);			
			unitCell.id_leak[i].blockId = ijk;
			unitCell.id_leak[i].order = (unitCell.co[ijk]);
		}

		/*for (int i = 0; i < unitCell.id_leak.size(); i++)
			unitCell.id_leak[i].order--;*/
		return unitCell;
	}
	else {
		std::cout << "Not ready:: Voro++ \n";
		return container_leak(1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1, 1, 1, 1);
	}
}

void VoronoiCellRoutine(const PeriodicCellList<Empty> &c, int numThreads, std::vector<double> * cellVolume, std::vector<int> * vertices, std::vector<int> * faces, std::vector<GeometryVector> * centroids) {
	bool VOL = (cellVolume != nullptr);
	bool VERT = (vertices != nullptr);
	bool FACE = (faces != nullptr);
	bool CENTR = (centroids != nullptr);
	PeriodicCellList<Empty> C(c);
	int d = C.GetDimension();
	int NUM = C.NumParticle();
	if (VOL) {
		(*cellVolume).clear();
		(*cellVolume).resize(NUM, 0.0);
	}
	if (VERT) {
		(*vertices).clear();
		(*vertices).resize(NUM, 0);
	}
	if (FACE) {
		(*faces).clear();
		(*faces).resize(NUM, 0);
	}

	if (CENTR) {
		(*centroids).clear();
		(*centroids).resize(NUM, GeometryVector(static_cast<int> (d)));
	}
	if (d > 1) {//USE VORO++
		if (VOL || VERT || FACE || CENTR) {
/*			int NUM = C.NumParticle();
			if (VOL) {
				(*cellVolume).clear();
				(*cellVolume).resize(NUM, 0.0);
			}
			if (VERT) {
				(*vertices).clear();
				(*vertices).resize(NUM, 0);
			}
			if (FACE) {
				(*faces).clear();
				(*faces).resize(NUM, 0);
			}

			if (CENTR) {
				(*centroids).clear();
				(*centroids).resize(NUM, GeometryVector(static_cast<int> (d)));
			}
*/
			int NMAX;
			double fac = 1;
			{
				switch (d) {
				case 2:
					NMAX = 10000;/*(int)round(100000);500 */  break;
				case 3:
					NMAX = (int)round(50000);/*10000000;*/ break;
				default:
					NMAX = 1;
					std::cout << "dimension is too high" << std::endl;
					break;
				}
			}
			int nx, ny, nz = 0;	bool trans;
			int mem = 30;
			int count_succ = 0, count_fail = 0;

			if (d == 2 || d == 3) {
				/*configuration is small enough to swallow it once*/
				if (NUM < NMAX) {
					std::vector<GeometryVector> basis;
					std::vector<container_leak> con_list;
					GetConversionParam(C, basis, trans, nx, ny, nz, mem);
					con_list.emplace_back(container_leak(C, basis, trans, nx, ny, nz, mem));
					/*Start computations!*/
					voro::voronoicell vc;
					voro::c_loop_all_periodic routine(con_list[0]);
					int idx = 0;

					std::vector<double> vrts;
					if (routine.start()) {
						do {
							idx = routine.pid();
							if (idx > -1) {/*Particles are in the inner region*/
								if (con_list[0].compute_cell(vc, routine)) {//Successfully computed
									//Cell Volume
									if (VOL) {
										(*cellVolume)[idx] = vc.volume();
									}
									//Cell Vertices
									if (VERT) {
										vc.vertices(vrts);
										int numV = vrts.size() / 3;
										if (d == 2)
											(*vertices)[idx] = numV / 2;
										else if (d == 3)
											(*vertices)[idx] = numV;
									}
									//Cell Faces
									if (FACE) {
										(*faces)[idx] = vc.number_of_faces();
									}
									//Cell centroids
									if (CENTR) {
										double x, y, z;
										vc.centroid(x, y, z);
										if (d == 2)
											(*centroids)[idx] = GeometryVector(x, y);
										else if (d == 3)
											(*centroids)[idx] = GeometryVector(x, y, z);
									}
									count_succ++;
								}
								else
									count_fail++;
							}
						} while (routine.inc());
					}
				}
				/*a configuration is too large => split the configuration*/
				else {
					double margin = 6.0;	double rho = C.NumParticle() / C.PeriodicVolume();
					mem = 10;
					int numSplit = (int)ceil(pow(C.PeriodicVolume(), 1.0 / d) / (pow(fac*NMAX / rho, 1.0 / d) - 2.0*margin));
					std::vector<GeometryVector> basis_;//Basis vectors for the domain of our interst in each split
					std::vector<GeometryVector> basis_subdomain;
					std::vector<GeometryVector> centerDomain; //In Cartesian coordinates.

					for (int i = 0; i < d; i++)
						basis_.push_back((1.0 / (double)numSplit) * C.GetBasisVector(i));
					{/*Center of domains*/
						GeometryVector center;
						for (int i = 0; i < numSplit; i++) {
							for (int j = 0; j < numSplit; j++) {
								if (d == 2) {
									center = (0.5 + i)*basis_[0] + (0.5 + j)*basis_[1];
									centerDomain.push_back(center);
								}
								else if (d == 3) {
									for (int k = 0; k < numSplit; k++) {
										center = (0.5 + i)*basis_[0] + (0.5 + j)*basis_[1] + (0.5 + k)*basis_[2];
										centerDomain.push_back(center);
									}
								}

							}
						}
					}
					for (int i = 0; i < d; i++) {//Extend the basis vectors of subdomains.
						basis_subdomain.push_back((1.0 + 2.0*margin / pow(rho, 1.0 / d) / sqrt(basis_[i].Modulus2()))*basis_[i]);
					}

					/*new basis vectors when basis_subdomain is not suitable for container class*/
					std::vector<GeometryVector> basis_subdomain_standard;
					std::vector<GeometryVector> unitbasisVectors;
					for (int i = 0; i < d; i++)
						unitbasisVectors.push_back((1.0 / sqrt(basis_[i].Modulus2()))*basis_[i]);

					Configuration subdomain_template(d, &basis_subdomain[0], 1.0);
					GetConversionParam(subdomain_template, (int)ceil(rho*subdomain_template.PeriodicVolume()), basis_subdomain_standard, trans, nx, ny, nz, mem);
					double r_max = 0.5 * subdomain_template.MaxDiagonalLength();
					C.PrepareIterateThroughNeighbors(1.01*r_max);

					//half lengths of basis vectors
					std::vector<double> region0;
					std::vector<double> region1;
					for (int i = 0; i < d; i++) {
						region0.push_back(0.501*sqrt(basis_[i].Modulus2()));
						region1.push_back(0.5*sqrt(basis_subdomain[i].Modulus2()));
					}
					/* return 0 = region of interest
					return 1 = marginal region
					return 2 = don't consider
					*/
					auto regionDet = [&unitbasisVectors, &region0, &region1](const GeometryVector & cartesianCoord) ->int {
						int region = 0;
						int region_idx;
						double projLength;
						for (int i = 0; i < unitbasisVectors.size(); i++) {
							projLength = fabs(unitbasisVectors[i].Dot(cartesianCoord));
							if (projLength <= region0[i])
								region_idx = 0;
							else if (projLength < region1[i])
								region_idx = 1;
							else
								region_idx = 2;

							region = std::max(region, region_idx);
						}

						return region;
					};
					GeometryVector ver2cen(basis_subdomain[0]);
					for (int i = 1; i < d; i++)
						ver2cen = ver2cen + basis_subdomain[i];
					ver2cen = 0.5*ver2cen;

					//Make a big (nonperiodic) rectangular box containing a single fundamental cell!
					double ax = 0, ay = 0, az = 0;	//lower bounds
					double bx = 0, by = 0, bz = 0;	//upper bounds
					ax = std::min(ax, basis_subdomain_standard[1].x[0]);
					ay = std::min(ay, basis_subdomain_standard[1].x[1]);
					bx = std::max(basis_subdomain_standard[0].x[0], basis_subdomain_standard[1].x[0]);
					by = std::max(basis_subdomain_standard[0].x[1], basis_subdomain_standard[1].x[1]);
					if (d == 2) {
						az = -0.5;
						bz = 0.5;
					}
					else if (d == 3) {
						ax = std::min(ax, basis_subdomain_standard[2].x[0]);
						ay = std::min(ay, basis_subdomain_standard[2].x[1]);
						bx = std::max(bx, basis_subdomain_standard[2].x[0]);
						by = std::max(by, basis_subdomain_standard[2].x[1]);
						az = std::min(az, basis_subdomain_standard[2].x[2]);
						bz = std::max(bz, basis_subdomain_standard[2].x[2]);
					}
					Configuration c_(d, &basis_subdomain_standard[0], 1);
#pragma omp parallel for  schedule(guided) num_threads(numThreads)
					for (int i = 0; i < centerDomain.size(); i++) {
						int numprts = 0;
						GeometryVector x0 = C.CartesianCoord2RelativeCoord(centerDomain[i]);
						//container_leak con(subdomain_template, basis_subdomain_standard, trans, nx, ny, nz, mem);
						voro::container con(ax, bx, ay, by, az, bz, nx, ny, nz, false, false, false, mem);
						//Input particles!
						C.IterateThroughNeighbors(x0, 1.01*r_max,
							[&con, &regionDet, &ver2cen, trans, &c_, &subdomain_template, &numprts](const GeometryVector & shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, size_t Sourceparticle)->void {
							int region = regionDet(shift);
							if (region <= 1) {
								GeometryVector r = shift + ver2cen;
								if (trans) {
									//r = con.c_.RelativeCoord2CartesianCoord(con.c0.CartesianCoord2RelativeCoord(r));
									r = c_.RelativeCoord2CartesianCoord(subdomain_template.CartesianCoord2RelativeCoord(r));
								}
								double x = r.x[0], y = r.x[1], z = 0.0;//z = 0.5;
								if (shift.Dimension == 3)
									z = r.x[2];
								int prt_idx = (region == 0) ? Sourceparticle : -1;
								con.put(prt_idx, x, y, z);	numprts++;
								//int ijk = con.put_ijk(prt_idx, x, y, z);
								//IDs temp;	temp.blockId = ijk;	temp.order = (con.co[ijk]) - 1;
								//con.id_leak.push_back(temp);
							}
						}
						);
						if (Verbosity > 5) {
#pragma omp critical
							{
								//std::cout << i << "\tnum prts = " << con.id_leak.size() << std::endl;
								std::cout << i << "\tnum prts = " << numprts << std::endl;
							}
						}
						//Compute Voronoi cells
						voro::voronoicell vc;
						//voro::c_loop_all_periodic routine(con);
						voro::c_loop_all routine(con);
						int idx = 0, count = 0, count_ = 0;
						std::vector<double> vrts;
						if (routine.start()) {
							do {
								idx = routine.pid();
								if (idx > -1) {/*Particles are in the inner region*/
									if (con.compute_cell(vc, routine)) {//Successfully computed
										if (VOL) {//Cell Volume
											(*cellVolume)[idx] = vc.volume();
										}
										if (VERT) {//Cell Vertices
											vc.vertices(vrts);
											int numV = vrts.size() / 3;
											if (d == 2)
												(*vertices)[idx] = numV / 2;
											else if (d == 3)
												(*vertices)[idx] = numV;
										}
										if (FACE) {//Cell Faces
											(*faces)[idx] = vc.number_of_faces();
										}
										if (CENTR) {//Cell centroids
											double x, y, z;
											vc.centroid(x, y, z);
											if (d == 2)
												(*centroids)[idx] = GeometryVector(x, y);
											else if (d == 3)
												(*centroids)[idx] = GeometryVector(x, y, z);
										}
										count++;
									}
									else
										count_++; //fail in calculation
								}
							} while (routine.inc());
						}
#pragma omp atomic
						count_fail += count_;
#pragma omp atomic
						count_succ += count;
					}
					if (Verbosity > 3) {
						std::cout << "# of success = " << count_succ << std::endl;
						std::cout << "# of failure = " << count_fail << std::endl;
					}
				}
			}
		}
		else {
			std::cout << "What do you want?\n";
		}
	}
	else if (d == 1) {
		if (VOL || VERT || FACE || CENTR) {
			//Align particle ID based on coordinates
			Merge_Sort(C, 0);
			//using the nearest neighbors on the right and left sides
			for (int i = 0; i < NUM; i++) {
				double x_beg=0.0, x_end=0.0;
				double x_cen = C.GetCartesianCoordinates(i).x[0];
				int idx_beg = i - 1, idx_end = i + 1;
				if (idx_beg < 0) {
					idx_beg = NUM - 1;
					x_beg = -C.GetBasisVector(0).x[0];
				}
				x_beg += C.GetCartesianCoordinates(idx_beg).x[0];

				if (idx_end == NUM) {
					idx_end = 0;
					x_end = +C.GetBasisVector(0).x[0];
				}
				x_end += C.GetCartesianCoordinates(idx_end).x[0];

				//find cell bounaries
				x_beg = 0.5*(x_beg + x_cen);
				x_end = 0.5*(x_end + x_cen);

				if (VOL)
					(*cellVolume)[i] = x_end - x_beg;
				if (CENTR)
					(*centroids)[i].x[0] = 0.5*(x_beg + x_end);
				if (VERT)
					(*vertices)[i] = 2;
				if (FACE)
					(*faces)[i] = 1;
			}

		}
		else
			std::cout << "What do you want?\n";
	}

}
void VoronoiCellRoutine_radial(const SpherePacking &c, int numThreads, std::vector<double> * cellVolume, std::vector<int> * vertices, std::vector<int> * faces, std::vector<GeometryVector> * centroids, bool radial) {
	if (!radial)//just centroidal voronoi tessellation
		VoronoiCellRoutine(c, numThreads, cellVolume, vertices, faces, centroids);
	else {
		//Please implement the radial Voronoi tessellation	
	}

}

#else
	void VoronoiCellRoutine(const PeriodicCellList<Empty> &c, std::vector<double> &cellVolume, std::vector<int> & vertices, std::vector<int> & faces, std::vector<GeometryVector> & centroids) {
		std::cout << "Please include \"VORO++\"\n";
	}
	void VoronoiCellRoutine_radial(const SpherePacking &c, std::vector<double> &cellVolume, std::vector<int> & vertices, std::vector<int> & faces, std::vector<GeometryVector> & centroids, bool radial = false) {
		std::cout << "Please include \"VORO++\"\n";
	}
#endif


#include <boost/math/special_functions/spherical_harmonic.hpp>
std::complex<double> Psi6(const PeriodicCellList<Empty> & Config, size_t j, int NeighborType, signed int m)
{
	typedef std::complex<double> complex;
	std::vector<GeometryVector> NeighborDisplacements;
	DimensionType dim = Config.GetDimension();
	if (NeighborType == 0)
	{
		std::vector<GeometryVector> neighbors;
		double TypicalLength = std::pow(Config.PeriodicVolume() / Config.NumParticle(), 1.0 / Config.GetDimension())*1.5;
		double l = TypicalLength;
		size_t Nneighbor;
		if (dim == 2)
			Nneighbor = 6;
		else if (dim == 3)
			Nneighbor = 12;
		else
			std::cerr << "Error in Psi6 : unsupported dimension!\n";
		while (neighbors.size() < Nneighbor+1)
		{
			neighbors.clear();
			Config.IterateThroughNeighbors(Config.GetRelativeCoordinates(j), l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
				neighbors.push_back(shift);
			});
			l *= 2;
		}
		auto CompareFunc = [](const GeometryVector & left, const GeometryVector & right) ->bool
		{
			return left.Modulus2() < right.Modulus2();
		};
		std::sort(neighbors.begin(), neighbors.end(), CompareFunc);
		for (int i = 1; i < Nneighbor + 1; i++)
			NeighborDisplacements.push_back(neighbors[i]);
	}
	else if (NeighborType == 1)
	{
		//std::vector<GeometryVector> neighbors;
		//DimensionType Dimension = Config.GetDimension();
		//std::vector<GeometryVector> vertices;
		//std::vector<signed long> Linkage;
		//GetVoronoiCell(Config, j, vertices, &Linkage);
		//for (size_t i = 0; i < vertices.size(); i++)
		//{
		//	for (DimensionType k = 0; k < Dimension + 1; k++)
		//	{
		//		if (Linkage[i*(Dimension + 1) + k] == -1)
		//			continue;
		//		if (Linkage[i*(Dimension + 1) + k] >= i)
		//			continue;
		//		GeometryVector cj = Config.GetCartesianCoordinates(j);
		//		GeometryVector vi = vertices[i];
		//		GeometryVector vk = vertices[Linkage[i*(Dimension + 1) + k]];
		//		GeometryVector cvi = vi - cj;
		//		GeometryVector vik = vk - vi;
		//		GeometryVector ortho = cvi - cvi.Dot(vik) / vik.Modulus2()*vik;
		//		neighbors.push_back(2.0*ortho);
		//	}
		//}

		//debug temp
		//std::cout << "Neighbors contain\n";
		//for (auto iter = neighbors.begin(); iter != neighbors.end(); ++iter)
		//	std::cout << *iter;

		//new code
		{
			DimensionType Dimension = Config.GetDimension();
			std::vector<GeometryVector> vertices;
			std::vector<signed long> Linkage;
			std::set<size_t> Neighbors;
			GetVoronoiCell(Config, j, vertices, &Linkage, &Neighbors);

			//debug temp
			//std::cout << "new Neighbors contain\n";
			for (auto iter = Neighbors.begin(); iter != Neighbors.end(); ++iter)
			{
				GeometryVector neighbor = Config.GetRelativeCoordinates(*iter) - Config.GetRelativeCoordinates(j);
				Config.RelativeCoordToMinimumImage(neighbor);
				neighbor = Config.RelativeCoord2CartesianCoord(neighbor);
				NeighborDisplacements.push_back(neighbor);
				//debug temp
				//std::cout << neighbor;
			}
		}

		//complex sum = 0;
		//for (int i = 0; i < neighbors.size(); i++)
		//	sum += std::exp(complex(0.0, 6.0)*std::atan2(neighbors[i].x[1], neighbors[i].x[0]));
		//return sum / (double)(neighbors.size());

	}
	else
	{
		std::cerr << "Error in Psi6 : unsupported NeighborType!\n";
		return complex(0, 0);
	}
	complex sum = 0;
	if (dim == 2)
		for (auto iter = NeighborDisplacements.begin(); iter != NeighborDisplacements.end(); ++iter)
			sum += std::exp(complex(0.0, 6.0)*std::atan2(iter->x[1], iter->x[0]));
	else if (dim == 3)
		for (auto iter = NeighborDisplacements.begin(); iter != NeighborDisplacements.end(); ++iter)
			sum += boost::math::spherical_harmonic(6, m, std::acos(iter->x[2] / std::sqrt(iter->Modulus2())), std::atan2(iter->x[1], iter->x[0]));
			//sum += std::exp(complex(0.0, m)*std::atan2(iter->x[1], iter->x[0])) * LegenderSphericalPlm(6, m, iter->x[2] / std::sqrt(iter->Modulus2()));
	else
		std::cerr << "Error in Psi6 : unsupported dimension!\n";
	return sum / (double)(NeighborDisplacements.size());
}

