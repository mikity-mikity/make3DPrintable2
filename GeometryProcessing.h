#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>

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

#include <vector>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include"face.h"
#include"halfedge.h"
#include"vertex.h"
#include<iostream>
using namespace std;
using namespace Eigen;
namespace GeometryProcessing
{
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
/*	struct Point3d{
	public:
		long double X,Y,Z;
		Point3d()
		{
			X=0;
			Y=0;
			Z=0;
		}
		Point3d(long double _x,long double _y,long double _z)
		{
			X=_x;
			Y=_y;
			Z=_z;
		}
		Point3d cross(Point3d f)
		{
			Point3d ret(Y*f.Z - Z*f.Y, Z*f.X - X*f.Z, X*f.Y - Y*f.X);
			return ret;
		}
		long double dot(Point3d f)
		{
			return X*f.X + Y*f.Y + Z*f.Z;
		}
		long double norm()
		{
			return std::sqrt(X*X + Y*Y + Z*Z);
		}
		Point3d operator-(Point3d d)
		{
			return Point3d(X - d.X, Y - d.Y, Z - d.Z);
		}
		Point3d operator/(long double f)
		{
			return Point3d(X / f, Y / f, Z / f);
		}
	};*/
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
		Eigen::SparseMatrix<halfedge*> __halfedgeTable;
		Eigen::SparseMatrix<vector<face*>*> _faceTable;
        orient* __orientation;
	public:
		vector<halfedge*> edges()
        {
            auto res = vector<halfedge*>();
            for (auto e : halfedges)
            {
                if (e->P->N < e->next->P->N)
                    res.push_back(e);
            }
            return res;
        }
        __int64 nVertices()
        {
            return vertices.size();
        }
        __int64 nFaces()
        {
            return faces.size();
        }
	private:
		MeshStructure()
        {
            this->Clear();
        }
        void Construct(Mesh *val)
        {
            int _nVertices = (int)val->Vertices.size();
            int _nFaces = (int)val->Faces.size();

            __orientation = new orient[_nFaces];
            _faceTable = SparseMatrix<vector<face*>*>(_nVertices, _nVertices);
            __halfedgeTable = SparseMatrix<halfedge*>(_nVertices, _nVertices);
			for (int i = 0; i < _nFaces; i++)
            {
                __orientation[i] = orient::unknown;
            }

            for (int i = 0; i < _nVertices; i++)
            {
                auto _v = new vertex(i);
                vertices.push_back(_v);
            }

            for (int i = 0; i < _nFaces; i++)
            {
                auto f = val->Faces[i];
                auto _f = new face(i, f.A, f.B, f.C);
                faces.push_back(_f);
                faceTableAdd(_f);
            }
            //Recursive
            halfEdgeAdd(faces[0]);
            //find pairs
            for (auto h : halfedges)
            {
                int i = h->P->N;
                int j = h->next->P->N;
				//if (__halfedgeTable.coeff(i, j) != NULL) throw new ArgumentOutOfRangeException(";)");
                __halfedgeTable.coeffRef(i, j) = h;
            }
            for (auto h : halfedges)
            {
                int i = h->P->N;
                int j = h->next->P->N;
                //if boundary edge...
                if (__halfedgeTable.coeff(j, i) == NULL)
                {
                    h->pair = NULL;
                }
                else
                {
                    h->pair = __halfedgeTable.coeffRef(j, i);
                }
            }
            //post process to find boundary vertices

			//align the first half edge at boundary
            for (auto v : vertices)
            {
                auto h = v->hf_begin;
                do
                {
                    if (h->prev->isNaked())
                    {
                        while (!h->isNaked())
                        {
                            h = h->pair->next;
                        }
                        v->hf_begin = h;
                        break;
                    }
                    h = h->prev->pair;
                } while (h != v->hf_begin);
            }

            vector<vector<halfedge*>> boundary_complements;
			int nBoundary = 0;
            for (auto v : vertices)
            {
                auto h = v->hf_begin;
				if (h->isNaked())//first naked halfedge found
				{
					bool flag = true;
					for (int i = 0; i < nBoundary; i++)
					{
						if (std::find(boundary_complements[i].begin(), boundary_complements[i].end(), h) != boundary_complements[i].end())
						{
							flag = false;
							break;
						}
					}
					if (flag)
					{
						boundary_complements.push_back(vector<halfedge*>());
						do
						{
							boundary_complements[nBoundary].push_back(h);
							h = h->next->P->hf_begin;
						} while (h != v->hf_begin);
						nBoundary++;
					}
				}
            }
            //insert boundary halfedges
			for (auto boundary_complement : boundary_complements)
			{
				vector<halfedge*> boundary;
				for (int i = 0; i < boundary_complement.size(); i++)
				{
					boundary.push_back(new halfedge(boundary_complement[i]->next->P));
					boundary[i]->pair = boundary_complement[i];
					boundary_complement[i]->pair = boundary[i];
					halfedges.push_back(boundary[i]);
				}
				boundaryStart.push_back(boundary[0]);
				for (int i = 0; i < boundary.size(); i++)
				{
					boundary[i]->owner = NULL;
					if (i != 0)
					{
						boundary[i]->next = boundary_complement[i - 1]->pair;
					}
					else
					{
						boundary[i]->next = boundary_complement[boundary_complement.size() - 1]->pair;
					}
					if (i != boundary.size() - 1)
					{
						boundary[i]->prev = boundary_complement[i + 1]->pair;
					}
					else
					{
						boundary[i]->prev = boundary_complement[0]->pair;
					}
				}
			}
            //check if any naked halfedge survives
            for (auto e : halfedges)
            {
                if (e->isNaked()) std::cout<<"naked halfedges survive!";
            }

			//post process to create stars
            for (auto v : vertices)
            {
                auto h = v->hf_begin;
                v->star.clear();
                do
                {
                    v->star.push_back(h);
                    if (h->isBoundary()) break;
                    h = h->prev->pair;
                } while (h != v->hf_begin);
            }
            //post process to create onering
            for (auto v : vertices)
            {
                auto h = v->hf_begin;
                v->onering.clear();
                do
                {
                    do
                    {
                        h = h->next;
                        v->onering.push_back(h);
                    } while (h->next->next->P != v);
                    if (h->next->pair->isBoundary()) break;
                    h = h->next->pair;
                } while (h != v->hf_begin);
            }
            //post process to split the vertices into inner and outer.
            innerVertices.clear();
            outerVertices.clear();
            for (auto v : vertices)
            {
                if (v->hf_begin->pair->isBoundary()) outerVertices.push_back(v); else innerVertices.push_back(v);
            }
			delete(__orientation);
        }
		void Construct_already_oriented(Mesh *val, vector<Delaunay::Facet> facet_list)
		{
			int _nVertices = (int)val->Vertices.size();
			int _nFaces = (int)val->Faces.size();
			faces.clear();
			vertices.clear();

			for (int i = 0; i < _nVertices; i++)
			{
				auto _v = new vertex(i);
				vertices.push_back(_v);
			}
			vector<halfedge*>* vectorPool = new vector<halfedge*>[val->Vertices.size()];
			int cc = 0;
			int tt = 0;
			for (int i = 0; i < _nFaces;i++)
			{
				auto f = val->Faces[i];
				auto _f = new face(i,f.A,f.B,f.C);
				faces.push_back(_f);
				halfedge* eA = new halfedge(vertices[f.A]);// &halfedgePool[i * 3 + 0];
				halfedge* eB = new halfedge(vertices[f.B]);
				halfedge* eC = new halfedge(vertices[f.C]);
				eA->next = eB;
				eB->next = eC;
				eC->next = eA;
				eA->owner = _f;
				eB->owner = _f;
				eC->owner = _f;
				eA->prev = eC;
				eB->prev = eA;
				eC->prev = eB;
				//pair is not determined here.
				halfedges.push_back(eA);
				halfedges.push_back(eB);
				halfedges.push_back(eC);
				vectorPool[f.A].push_back(eA);
				vectorPool[f.B].push_back(eB);
				vectorPool[f.C].push_back(eC);
				//if ((i / 10000) * 10000 == i)std::cout << i << endl;
			}
			//circulate vertices and count up errors
			//find pairs
			int count1 = 0;
			int count2 = 0;
			int CCount1 = 0;
			int CCount2 = 0;
			double PI2 = boost::math::constants::pi<double>() * 2.;
			double PI = boost::math::constants::pi<double>();
			vector<int> errors;

			for (int i = 0; i < _nVertices; i++)
			{
				auto pool = vectorPool[i];
				vertex* v = vertices[i];
				bool flag = false;
				for (auto p : pool) //p is a halfedge
				{
					vertex* w=p->next->P;
					vector<std::pair<halfedge*, bool>> candidates;
					for (auto t : vectorPool[w->N])
					{
						if (t->next->P->N == i)
							candidates.push_back(std::make_pair(t,false));
					}
					/*for (auto t : pool)
					{
						if (t != p)
						{
							if (t->P->N == i)
							{
								if (t->next->P==w)
									candidates.push_back(std::make_pair(t, true));
							}
						}
					}*/
					if (candidates.size() == 0){
						std::cout << "candidates.size()==0" << endl;
						std::cin.get();
						exit(1);
					}else if (candidates.size()>1){
						count2++;
						face* FA = p->owner;
						vector<std::pair<halfedge*,bool>> candidates2;
						for (auto itr = candidates.begin(); itr != candidates.end(); itr++)
						{
							auto hf = itr->first;
							face* FB = hf->owner;
							Vector NA = val->Faces[FA->N].N;
							Vector NB = val->Faces[FB->N].N;
							Vector a = val->Vertices[p->P->N] - val->Vertices[p->prev->P->N];
							Vector d = val->Vertices[hf->next->next->P->N] - val->Vertices[hf->next->P->N];
							Vector V = val->Vertices[p->next->P->N] - val->Vertices[p->P->N];
							if (itr->second){
								d = val->Vertices[hf->prev->P->N] - val->Vertices[hf->P->N];
							}
							a = a - V*(a*V / V.squared_length());
							d = d - V*(d*V / V.squared_length());
							double valB = CGAL::cross_product(a, NA)*CGAL::cross_product(d, NB);
							if (valB > 0){
								candidates2.push_back(*itr);
							}
						}

						if (candidates2.size() == 1)
						{
							p->pair = candidates2[0].first;
							if (candidates2[0].second)
							{
								CCount2++;
							}
							else CCount1++;
						}
						else if (candidates2.size() == 0)
						{
							std::cout << "found an error" << endl;
							std::cout << "Press enter to continue" << endl;
							std::cin.get();
						}else
						{
							//calculate dihedral angle
							vector<double> angles;
							for (auto itr = candidates2.begin(); itr != candidates2.end(); itr++)
							{
								auto hf = itr->first;
								face* FB = hf->owner;
								Vector NA=val->Faces[FA->N].N;
								Vector NB=val->Faces[FB->N].N;
								Vector a = val->Vertices[p->P->N] - val->Vertices[p->prev->P->N];
								Vector d = val->Vertices[hf->next->next->P->N] - val->Vertices[hf->next->P->N];
								Vector V = val->Vertices[p->next->P->N] - val->Vertices[p->P->N];
								if (itr->second){
									d = val->Vertices[hf->prev->P->N] - val->Vertices[hf->P->N];
								}
								a = a - V*(a*V / V.squared_length());
								d = d - V*(d*V / V.squared_length());

								double val = CGAL::cross_product(a, NA)*CGAL::cross_product(a, d);
								bool leftorright = false;
								if (val < 0)leftorright = true;  //if true, it is turning to the right...
								//calculate angle
								double cos = NA*NB;
								double sin = std::sqrt(CGAL::cross_product(NA, NB).squared_length());
								if (cos > 1.)cos = 1.;
								if (cos < -1.)cos = -1.;
								if (sin > 1.)sin = 1.;
								if (sin < -1.)sin = -1.;
								if (!leftorright)sin = -sin;
								double theta1 = std::acos(cos);  //0-180
								double theta2 = std::asin(sin);  //-90-90
								double theta = 0;
								if (sin >= 0 && cos >= 0)theta = theta1;//0-90
								if (sin >= 0 && cos < 0)theta = theta1;//90-180
								if (sin < 0 && cos < 0)theta = -theta1;//-180 - -90
								if (sin < 0 && cos >= 0)theta = theta2;//-90 - 0
								angles.push_back(theta);
							}
							//choose index of max theta
							double max = -1000;
							int maxIndex = -1000;
							for (int k = 0; k < angles.size(); k++)
							{
								if (angles[k]>max){
									max = angles[k];
									maxIndex = k;
								}
							}
							p->pair = candidates2[maxIndex].first;
							if (candidates2[maxIndex].second)
							{
								CCount2++;
							}
							else CCount1++;
						}
					}
					else
					{
						p->pair = candidates[0].first;
						count1++;
						if (candidates[0].second)
						{
							CCount2++;
						}
						else CCount1++;
					}
				}
			}
			std::cout << "count1:" << count1 << "/" << "count2:" << count2 << endl;
			std::cout << "CCount1:" << CCount1 << "/" << "CCount2:" << CCount2 << endl;
			count1 = 0;
			count2 = 0;
			int count3 = 0;
			for (auto hf : halfedges)
			{
				if (hf->pair != NULL)
				{
					if (hf->pair->pair != NULL)
					{

						if (hf->pair->pair == hf)
						{
							count1++;
						}
						else
						{
							count2++;
							std::cout << "error:" << hf->P->N << "-" << hf->next->P->N << endl;
						}
					}
					else
					{
						std::cout << "hf->pair->pair=NULL" << endl;
						count3++;
					}
				}
				else
				{
					std::cout << "hf->pair=NULL" << endl;
					count3++;
				}
			}
			std::cout << "count1:" << count1 << "/" << "count2:" << count2 << "/" << "count3:" << count3 << endl;
			//onestars
			count1 = 0;
			count2 = 0;
			count3 = 0;
			int __nVertices = _nVertices;
			for (int i = 0; i < __nVertices; i++)
			{
				auto pool = vectorPool[i];
				auto v = vertices[i];
				vector<vector<halfedge*>>stars;
				do
				{
					auto h = pool[0];
					vector<halfedge*>star;
					do{
						star.push_back(h);
						pool.erase(std::find(pool.begin(), pool.end(), h));
						if (h->isBoundary()) break;
						h = h->prev->pair;
					} while (h != star[0]);
					stars.push_back(star);
				} while (!pool.empty());
				if (stars.size() == 1){
					count1++;
					if (stars[0].size() == 0)
					{
						std::cout << "error" << endl;
						std::cin.get();
					}
					v->hf_begin = stars[0][0];
					v->star = stars[0];
				}
				double x, y, z;
				if (stars.size() >= 2){
					count2++;
					v->hf_begin = stars[0][0];
					v->star = stars[0];
					Vector VV(0, 0, 0);

					for (auto s : v->star)
					{
						auto cell=facet_list[s->owner->N].first;
						auto D = facet_list[s->owner->N].second;
						auto V=cell->vertex(D)->point()-val->Vertices[i];
						VV = VV+V;
					}
					VV = VV / v->star.size();
					VV = VV/1000;
					val->Vertices[i] = val->Vertices[i] + VV;
					int SS = stars.size() - 1;
					for (int j = 0; j < SS; j++)
					{
						vertices.push_back(new vertex(_nVertices + j));
						VV=Vector(0, 0, 0);
						for (auto s : stars[j + 1])
						{
							auto cell = facet_list[s->owner->N].first;
							auto D = facet_list[s->owner->N].second;
							auto V = cell->vertex(D)->point() - val->Vertices[i];
							VV = VV + V;
						}
						VV = VV / v->star.size();
						VV = VV / 1000;
						val->Vertices.push_back(val->Vertices[i] + VV);
						vertices[_nVertices + j]->hf_begin = stars[j + 1][0];
						vertices[_nVertices + j]->star = stars[j + 1];

						for (auto h : stars[j + 1])
						{
							h->P = vertices[_nVertices + j];
							for (int k = 0; k < 3; k++)
							{
								if (h->owner->corner[k] == v->N){
									h->owner->corner[k] = _nVertices + j;
									if (k==0)
										val->Faces[h->owner->N].A = _nVertices + j;
									if (k == 1)
										val->Faces[h->owner->N].B = _nVertices + j;
									if (k == 2)
										val->Faces[h->owner->N].C = _nVertices + j;
								}
							}
						}
					}
					_nVertices += SS;
				}
				if (i == 113611)std::cin.get();

			}
			std::cout << "stars" << endl;
			std::cout << "count1:" << count1 << "/" << "count2:" << count2 << endl;
			//post process to create onering
			for (auto v : vertices)
			{
				auto h = v->hf_begin;
				v->onering.clear();
				do
				{
					do
					{
						h = h->next;
						v->onering.push_back(h);
					} while (h->next->next->P != v);
					if (h->next->pair->isBoundary()) break;
					h = h->next->pair;
				} while (h != v->hf_begin);
			}
			/*			for (int i = 0; i < val->Vertices.size(); i++)
			{
				vector<vector<halfedge*>> onestars;
				auto pool = vectorPool[i];
				vertex* v = vertices[i];
				while (pool.size() != 0)
				{
					vector<halfedge*> onestar;
					onestar.push_back(pool.back());
					pool.pop_back();
					onestars.push_back(onestar);
				}
			}*/
/*			int count = 0;
			vector<int> errors;
			for(int i=0;i<val->Vertices.size();i++)
			{
				auto pool=vectorPool[i];
				vector<int> mm;
				for(auto p:pool)
				{
					mm.push_back(p->next->P->N);
				}
				bool flag = true;
				
				std::set<int> set;
				//std::vector<int> result;
				for (auto it = mm.begin(); it !=mm.end(); ++it)
				{
					if (!set.insert(*it).second)
						flag = false;//result.push_back(*it);
				}
				if (!flag){
					count++;
					errors.push_back(i);
				}
			}
			std::cout << count << "/" << val->Vertices.size() << endl;
			for (int i : errors)
			{
				auto pool = vectorPool[i];
				auto v = vertices[i];
				//fid pairs
			}*/
			//find pairs
			/*for (auto h : halfedges)
			{
				int i = h->P->N;
				int j = h->next->P->N;
				//if (__halfedgeTable.coeff(i, j) != NULL) throw new ArgumentOutOfRangeException(";)");
				__halfedgeTable.coeffRef(i, j) = h;
			}*/
			/*for (auto h : halfedges)
			{
				int i = h->P->N;
				int j = h->next->P->N;
				//if boundary edge...
				if (__halfedgeTable.coeff(j, i) == NULL)
				{
					h->pair = NULL;
				}
				else
				{
					h->pair = __halfedgeTable.coeffRef(j, i);
				}
			}*/
			//post process to find boundary vertices

			//align the first half edge at boundary
			/*for (auto v : vertices)
			{
				auto h = v->hf_begin;
				do
				{
					if (h->prev->isNaked())
					{
						while (!h->isNaked())
						{
							h = h->pair->next;
						}
						v->hf_begin = h;
						break;
					}
					h = h->prev->pair;
				} while (h != v->hf_begin);
			}

			vector<vector<halfedge*>> boundary_complements;
			int nBoundary = 0;
			for (auto v : vertices)
			{
				auto h = v->hf_begin;
				if (h->isNaked())//first naked halfedge found
				{
					bool flag = true;
					for (int i = 0; i < nBoundary; i++)
					{
						if (std::find(boundary_complements[i].begin(), boundary_complements[i].end(), h) != boundary_complements[i].end())
						{
							flag = false;
							break;
						}
					}
					if (flag)
					{
						boundary_complements.push_back(vector<halfedge*>());
						do
						{
							boundary_complements[nBoundary].push_back(h);
							h = h->next->P->hf_begin;
						} while (h != v->hf_begin);
						nBoundary++;
					}
				}
			}
			//insert boundary halfedges
			for (auto boundary_complement : boundary_complements)
			{
				vector<halfedge*> boundary;
				for (int i = 0; i < boundary_complement.size(); i++)
				{
					boundary.push_back(new halfedge(boundary_complement[i]->next->P));
					boundary[i]->pair = boundary_complement[i];
					boundary_complement[i]->pair = boundary[i];
					halfedges.push_back(boundary[i]);
				}
				boundaryStart.push_back(boundary[0]);
				for (int i = 0; i < boundary.size(); i++)
				{
					boundary[i]->owner = NULL;
					if (i != 0)
					{
						boundary[i]->next = boundary_complement[i - 1]->pair;
					}
					else
					{
						boundary[i]->next = boundary_complement[boundary_complement.size() - 1]->pair;
					}
					if (i != boundary.size() - 1)
					{
						boundary[i]->prev = boundary_complement[i + 1]->pair;
					}
					else
					{
						boundary[i]->prev = boundary_complement[0]->pair;
					}
				}
			}
			//check if any naked halfedge survives
			for (auto e : halfedges)
			{
				if (e->isNaked()) std::cout << "naked halfedges survive!";
			}

			//post process to split the vertices into inner and outer.
			innerVertices.clear();
			outerVertices.clear();
			for (auto v : vertices)
			{
				if (v->hf_begin->pair->isBoundary()) outerVertices.push_back(v); else innerVertices.push_back(v);
			}
			delete(__orientation);
			*/
		}
		private:
			void halfEdgeAdd(face *f)
			{
				auto _o = orient::unknown;
				for (int i = 0; i < 3; i++)
				{
					__int64 I = f->corner[i];
					__int64 J = 0;
					if(i == 2) J=f->corner[0];else J=f->corner[i + 1];
					if (_faceTable.coeffRef(I, J)->size() == 2)
					{
						if (_faceTable.coeffRef(I, J)->at(0) == f)
						{
							if (__orientation[_faceTable.coeffRef(I, J)->at(1)->N] != orient::unknown)
							{
								_o = orient::unknown;
								if(__orientation[_faceTable.coeffRef(I, J)->at(1)->N] == orient::positive)
									_o=orient::negative;else orient::positive;
							}
						}
						if (_faceTable.coeffRef(I, J)->at(1) == f)
						{
							if (__orientation[_faceTable.coeffRef(I, J)->at(0)->N] != orient::unknown)
							{
								_o = orient::unknown;
								if(__orientation[_faceTable.coeffRef(I, J)->at(0)->N] == orient::positive)
									_o=orient::negative;else _o=orient::positive;
							}
						}
					}
					else
					{
						if (_faceTable.coeff(J, I) != NULL)
						{
							if (__orientation[_faceTable.coeffRef(J, I)->at(0)->N] != orient::unknown)
							{
								_o = __orientation[_faceTable.coeffRef(J, I)->at(0)->N];
							}
						}
					}
				}
				__orientation[f->N] = orient::unknown;
				if(_o == orient::unknown)
					__orientation[f->N]=orient::positive;else __orientation[f->N]=_o;
				//register a halfedge
				if (__orientation[f->N] == orient::positive)
				{
					auto he1 = new halfedge(vertices[f->corner[0]]);
					auto he2 = new halfedge(vertices[f->corner[1]]);
					auto he3 = new halfedge(vertices[f->corner[2]]);
					halfedges.push_back(he1);
					halfedges.push_back(he2);
					halfedges.push_back(he3);
					he1->prev = he3; he1->next = he2; he1->owner = f;
					he2->prev = he1; he2->next = he3; he2->owner = f;
					he3->prev = he2; he3->next = he1; he3->owner = f;
					f->firsthalfedge = he1;
				}

				if (__orientation[f->N] == orient::negative)
				{
					auto he1 = new halfedge(vertices[f->corner[2]]);
					auto he2 = new halfedge(vertices[f->corner[1]]);
					auto he3 = new halfedge(vertices[f->corner[0]]);
					halfedges.push_back(he1);
					halfedges.push_back(he2);
					halfedges.push_back(he3);
					he1->prev = he3; he1->next = he2; he1->owner = f;
					he2->prev = he1; he2->next = he3; he2->owner = f;
					he3->prev = he2; he3->next = he1; he3->owner = f;
					f->firsthalfedge = he1;
				}


				//list up neighbors that are not oriented
				for (int i = 0; i < 3; i++)
				{
					__int64 I = f->corner[i];
					__int64 J = 0;
					if(i == 2) J=f->corner[0];else J=f->corner[i + 1];
					if (_faceTable.coeffRef(I, J)->size() == 2)
					{
						if (_faceTable.coeffRef(I, J)->at(0) == f)
						{
							if (__orientation[_faceTable.coeffRef(I, J)->at(1)->N] == orient::unknown)
							{
								halfEdgeAdd(_faceTable.coeffRef(I, J)->at(1));
							}
						}
						if (_faceTable.coeffRef(I, J)->at(1) == f)
						{
							if (__orientation[_faceTable.coeffRef(I, J)->at(0)->N] == orient::unknown)
							{
								halfEdgeAdd(_faceTable.coeffRef(I, J)->at(0));
							}
						}
					}
					else
					{
						if (_faceTable.coeff(J, I) != NULL)
						{
							if (__orientation[_faceTable.coeffRef(J, I)->at(0)->N] == orient::unknown)
							{
								halfEdgeAdd(_faceTable.coeffRef(J, I)->at(0));
							}
						}
					}
				}
			}
		private:
			void faceTableAdd(int i, int j, face* f)
			{
				if (_faceTable.coeff(i, j)==NULL)
				{
					_faceTable.coeffRef(i, j) = new vector<face*>();
				}
				_faceTable.coeffRef(i, j)->push_back(f);
			}
			void faceTableAdd(face* f)
			{
				for (int i = 0; i < 3; i++)
				{
					__int64 I = f->corner[i];
					__int64 J=0;
					if(i == 2) J=f->corner[0];else J= f->corner[i + 1];
					faceTableAdd(I, J, f);
				}
			}
		public:
			static MeshStructure* CreateFrom(Mesh *val)
			{
				MeshStructure* ret=new MeshStructure();
				ret->Construct(val);
				return ret;
			}
			static MeshStructure* CreateFrom_already_oriented(Mesh *val, vector<Delaunay::Facet> facet_list)
			{
				MeshStructure* ret = new MeshStructure();
				ret->Construct_already_oriented(val,facet_list);
				return ret;
			}
		public:
			void Clear()
			{
				if(vertices.size()!=0)
				{
					for(auto v : vertices)
					{
						delete(v);
					}
				}
				vertices.clear();
				if(faces.size()!=0)
				{
					for(auto f:faces)
					{
						delete(f);
					}
				}
				faces.clear();
				if(halfedges.size()!=0)
				{
					for(auto hf:halfedges)
					{
						delete(hf);
					}
				}
				if(_faceTable.nonZeros()!=0)
				{
					for (int k=0; k<_faceTable.outerSize(); ++k)
					{
					  for (SparseMatrix<vector<face*>*>::InnerIterator it(_faceTable,k); it; ++it)
					  {
						delete(it.value());						
						//it.row();   // row index
						//it.col();   // col index (here it is equal to k)
						//it.index(); // inner index, here it is equal to it.row()
					  }
					}
				}
				_faceTable.setZero();
				__halfedgeTable.setZero();
				halfedges.clear();
				innerVertices.clear();
				outerVertices.clear();
				boundaryStart.clear();
			}
			~MeshStructure()
			{
				this->Clear();
			}
    };
}