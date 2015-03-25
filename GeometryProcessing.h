#ifndef _GEOM_PROC_H_
#define _GEOM_PROC_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include<CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>


#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <CGAL/property_map.h>
#
// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

#include <vector>
//#include<Eigen/Dense>
//#include<Eigen/Sparse>
#include"face.h"
#include"halfedge.h"
#include"vertex.h"
#include<iostream>

namespace SMS = CGAL::Surface_mesh_simplification;

typedef CGAL::Simple_cartesian<long double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Poly;
typedef Poly::Point_3 PP;
typedef Poly::HalfedgeDS HalfedgeDS;

typedef boost::graph_traits<Poly>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Poly>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Poly>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<Poly>::edge_iterator edge_iterator;
struct Constrained_edge_map : public boost::put_get_helper<bool,Constrained_edge_map>
{
  typedef boost::readable_property_map_tag      category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef edge_descriptor                       key_type;

  Constrained_edge_map(const CGAL::Unique_hash_map<key_type,bool>& aConstraints)
    : mConstraints(aConstraints)
  {}

  reference operator[](key_type const& e) const { return  is_constrained(e); }

  bool is_constrained( key_type const& e ) const {
    return mConstraints.is_defined(e);
  }

private:
  const CGAL::Unique_hash_map<key_type,bool>& mConstraints;
};
struct col_int
{
	CGAL::Color col;
	__int64 num;
	col_int()
	{
		col = CGAL::BLACK;
		num = 0;
	}
	col_int(CGAL::Color _col, __int64 _num)
	{
		col = _col;
		num = _num;
	}


};
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<col_int, K> Vb;
typedef CGAL::Triangulation_cell_base_with_info_3<__int64, K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                      Delaunay;
//typedef CGAL::Delaunay_triangulation_3<K,CGAL::Fast_location> Delaunay;
typedef Delaunay::Point Point;
typedef CGAL::Vector_3<K> Vector;

using namespace std;
using namespace Eigen;
namespace GeometryProcessing
{
	bool is_border(edge_descriptor e, const Poly& sm);

	PP point(vertex_descriptor vd, const Poly& sm);


	struct MeshFace
	{
	public:
		__int64 A,B,C;
		Vector N;
		MeshFace(__int64 _A, __int64 _B, __int64 _C)
		{
			A = _A;
			B = _B;
			C = _C;
			N = Vector(0, 0, 0);
		}
		MeshFace(__int64 _A, __int64 _B, __int64 _C,Vector _N)
		{
			A = _A;
			B = _B;
			C = _C;
			N = _N;
		}
	};

	struct Mesh
	{
	public:
		vector<MeshFace> Faces;
		vector<Point> Vertices;
		Mesh()
		{
			Faces.clear();
			Vertices.clear();
		}
	};


    class MeshStructure
    {
	private:
		enum orient
        {
            unknown, positive, negative
        };

        //to get boundary chain
        //boundaryStart->hf->next->P->hf->next->P->....
	public:
        vector<halfedge*> boundaryStart;
        vector<vertex*> vertices;
        vector<face*> faces;
        vector<halfedge*> halfedges;
        vector<vertex*> innerVertices;
        vector<vertex*> outerVertices;  
	private:
		//halfedge** __halfedgeTable;
		//Eigen::SparseMatrix<halfedge*> __halfedgeTable;
		//Eigen::SparseMatrix<vector<face*>*> _faceTable;
		vector<halfedge*>* __halfedgeTable;
		vector<std::pair<face*, __int64>>* _faceTable;
		orient* __orientation;
	public:
		vector<halfedge*> edges();
		__int64 nVertices();
		__int64 nFaces();
	private:
		MeshStructure();
		void Construct(Mesh *val);
		void Construct_already_oriented(Mesh *val, vector<Delaunay::Facet> facet_list);
		void halfEdgeAdd(face *f);
		void faceTableAdd(int i, int j, face* f);
		void faceTableAdd(face* f);
	public:
		static MeshStructure* CreateFrom(Mesh *val);
		static Poly CreateFrom_already_oriented(Mesh *val, vector<Delaunay::Facet> facet_list);
		void Clear();
		~MeshStructure();
    };
	template<class HDS> class Builder:public CGAL::Modifier_base < HDS >{
	public:
		Mesh* mesh;
		MeshStructure *MS;
		void operator()(HDS& hds);
	};

}
#endif