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
		int A,B,C;
		double nx, ny, nz;
		MeshFace(int _A, int _B, int _C)
		{
			A = _A;
			B = _B;
			C = _C;
			nx = 0;
			ny = 0;
			nz = 0;
		}
		MeshFace(int _A, int _B, int _C,double _nx,double _ny,double _nz)
		{
			A = _A;
			B = _B;
			C = _C;
			nx = _nx;
			ny = _ny;
			nz = _nz;
		}
	};
	struct Point3d{
	public:
		double X,Y,Z;
		Point3d()
		{
			X=0;
			Y=0;
			Z=0;
		}
		Point3d(double _x,double _y,double _z)
		{
			X=_x;
			Y=_y;
			Z=_z;
		}
	};

	struct Mesh
	{
	public:
		vector<MeshFace> Faces;
		vector<Point3d> Vertices;
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
		void Construct_already_oriented(Mesh *val)
		{
			int _nVertices = (int)val->Vertices.size();
			int _nFaces = (int)val->Faces.size();


			for (int i = 0; i < _nVertices; i++)
			{
				auto _v = new vertex(i);
				vertices.push_back(_v);
			}
			int i = 0;
			vector<halfedge*>* vectorPool = new vector<halfedge*>[val->Vertices.size()];
			int cc = 0;
			for (auto f : val->Faces)
			{
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
				i++;
				if ((i / 10000) * 10000 == i)std::cout << i << endl;
			}			//find pairs
			/*for (auto h : halfedges)
			{
				int i = h->P->N;
				int j = h->next->P->N;
				//if (__halfedgeTable.coeff(i, j) != NULL) throw new ArgumentOutOfRangeException(";)");
				__halfedgeTable.coeffRef(i, j) = h;
			}*/
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
				if (e->isNaked()) std::cout << "naked halfedges survive!";
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
		private:
			void halfEdgeAdd(face *f)
			{
				auto _o = orient::unknown;
				for (int i = 0; i < 3; i++)
				{
					int I = f->corner[i];
					int J = 0;
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
					int I = f->corner[i];
					int J = 0;
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
					int I = f->corner[i];
					int J=0;
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
			static MeshStructure* CreateFrom_already_oriented(Mesh *val)
			{
				MeshStructure* ret = new MeshStructure();
				ret->Construct_already_oriented(val);
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