#include"GeometryProcessing.h"
#include"vertex.h"
#include"halfedge.h"
#include"face.h"
#include <boost/math/constants/constants.hpp>

namespace GeometryProcessing
{
	template<class HDS> void Builder<HDS>::operator()(HDS& hds){
		//Mesh* mesh;
		//MeshStructure *MS;
		CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
		B.begin_surface(mesh->Vertices.size(), mesh->Faces.size(), MS->halfedges.size());
		typedef typename HDS::Vertex Vertex;
		typedef typename Vertex::Point PPP;
		for (auto v : mesh->Vertices)
		{
			B.add_vertex(PPP(v.x(), v.y(), v.z()));
		}
		for (auto f : mesh->Faces)
		{
			B.begin_facet();
			B.add_vertex_to_facet(f.A);
			B.add_vertex_to_facet(f.B);
			B.add_vertex_to_facet(f.C);
			B.end_facet();
		}
		B.end_surface();
	}
	vector<halfedge*> MeshStructure::edges()
	{
		auto res = vector<halfedge*>();
		for (auto e : halfedges)
		{
			if (e->P->N < e->next->P->N)
				res.push_back(e);
		}
		return res;
	}
	__int64 MeshStructure::nVertices()
	{
		return vertices.size();
	}
	__int64 MeshStructure::nFaces()
	{
		return faces.size();
	}
	MeshStructure::MeshStructure()
	{
		this->Clear();
	}
	void MeshStructure::getfaces(int I, int J, vector<face*> &faces)
	{
		faces.clear();
		auto ff = _faceTable[I];
		for(auto f : ff)
		{
			if (f.second == J) faces.push_back(f.first);
		}
	}
	void MeshStructure::gethalfedges(int I, int J, vector<halfedge*> &lhalfedges)
	{
		lhalfedges.clear();
		auto ff = __halfedgeTable[I];
		for(auto f : ff)
		{
			if (f->next->P->N == J) lhalfedges.push_back(f);
		}
	}
	void MeshStructure::Construct(Mesh *val)
	{
		int _nVertices = (int)val->Vertices.size();
		int _nFaces = (int)val->Faces.size();

		__orientation = new orient[_nFaces];
		_faceTable = new vector<std::pair<face*,__int64>>[_nVertices];
		__halfedgeTable = new vector<halfedge*>[_nVertices];
		//_faceTable = SparseMatrix<vector<face*>*>(_nVertices, _nVertices);
		//__halfedgeTable = SparseMatrix<halfedge*>(_nVertices, _nVertices);
		for (int i = 0; i < _nVertices; i++)
		{
			_faceTable[i].clear();// = new List<tuple>();
			__halfedgeTable[i].clear();// = new List<halfedge>();
		}
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
		while (!neighbors.empty())
		{
			halfEdgeAdd(neighbors);
		}
		vector<halfedge*> lhalfedges;
		//find pairs
		for (auto h : halfedges)
		{
			int i = h->P->N;
			int j = h->next->P->N;
			gethalfedges(i, j, lhalfedges);
			if (lhalfedges.size() != 0) std::cout << "this should not be empty!" << endl;

			//if (__halfedgeTable.coeff(i, j) != NULL) throw new ArgumentOutOfRangeException(";)");
			//__halfedgeTable.coeffRef(i, j) = h;
			__halfedgeTable[i].push_back(h);

		}
		for (auto h : halfedges)
		{
			int i = h->P->N;
			int j = h->next->P->N;
			gethalfedges(j, i, lhalfedges);
			//if boundary edge...
			if (lhalfedges.size() == 0)
			{
				h->pair = NULL;
			}
			else
			{
				h->pair = lhalfedges[0];
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
		delete[](__orientation);
		delete[](_faceTable);
		delete[](__halfedgeTable);

	}
	void MeshStructure::Construct_already_oriented(Mesh *val, vector<Delaunay::Facet> facet_list)
	{
		__int64 _nVertices = (int)val->Vertices.size();
		__int64 _nFaces = (int)val->Faces.size();
		faces.clear();
		vertices.clear();

		for (int i = 0; i < _nVertices; i++)
		{
			auto _v = new vertex(i);
			vertices.push_back(_v);
		}
		vector<vector<halfedge*>> vectorPool;// = new vector<halfedge*>[val->Vertices.size()];
		for (int i = 0; i < val->Vertices.size(); i++)
		{
			vectorPool.push_back(vector<halfedge*>());
		}

		int cc = 0;
		int tt = 0;
		for (int i = 0; i < _nFaces; i++)
		{
			auto f = val->Faces[i];
			auto _f = new face(i, f.A, f.B, f.C);
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
				vertex* w = p->next->P;
				vector<std::pair<halfedge*, bool>> candidates;
				for (auto t : vectorPool[w->N])
				{
					if (t->next->P->N == i)
						candidates.push_back(std::make_pair(t, false));
				}
				for (auto t : pool)
				{
					if (t != p)
					{
						if (t->P->N == i)
						{
							if (t->next->P == w)
								candidates.push_back(std::make_pair(t, true));
						}
					}
				}
				if (candidates.size() == 0){
					std::cout << "candidates.size()==0" << endl;
					std::cin.get();
					exit(1);
				}
				else if (candidates.size() > 1){
					count2++;
					face* FA = p->owner;
					vector<std::pair<halfedge*, bool>> candidates2;
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
					}
					else
					{
						//calculate dihedral angle
						vector<double> angles;
						for (auto itr = candidates2.begin(); itr != candidates2.end(); itr++)
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
						if ((p->P->N == 129958 && p->next->P->N == 129959) || (p->P->N == 129959 && p->next->P->N == 129958))
						{
							std::cout << p->P->N << "," << p->next->P->N << "," << p->next->next->P->N << ":" << p->owner->N << endl;
							std::cout << "angle:"<<max<<endl;
							std::cout << p->pair->P->N << "," << p->pair->next->P->N << "," << p->pair->next->next->P->N << ":" << p->pair->owner->N << endl;
						}
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
		int __nVertices = _nVertices;
		
		//std::cin.get();

		


		//onestars
		count1 = 0;
		count2 = 0;
		count3 = 0;
		__nVertices = _nVertices;

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
			if (stars.size() >= 2){
				count2++;
				v->hf_begin = stars[0][0];
				v->star = stars[0];
				Vector VV(0, 0, 0);

				for (auto s : v->star)
				{
					auto cell = facet_list[s->owner->N].first;
					auto D = facet_list[s->owner->N].second;
					auto V = cell->vertex(D)->point() - val->Vertices[i];
					VV = VV + V;
				}
				VV = VV / ((double)v->star.size());
				VV = VV / ((double)10.);
				val->Vertices[i] = val->Vertices[i] + VV;
				__int64 SS = stars.size() - 1;
				for (__int64 j = 0; j < SS; j++)
				{
					vertices.push_back(new vertex(_nVertices + j));
					VV = Vector(0, 0, 0);
					for (auto s : stars[j + 1])
					{
						auto cell = facet_list[s->owner->N].first;
						auto D = facet_list[s->owner->N].second;
						auto V = cell->vertex(D)->point() - val->Vertices[i];
						VV = VV + V;
					}
					VV = VV / ((double)v->star.size());
					VV = VV / ((double)10.);
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
								if (k == 0)
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
		}
		
		std::cout << "stars" << endl;
		std::cout << "count1:" << count1 << "/" << "count2:" << count2 << "count3:" << count3 << endl;

		std::cout << "check shared halfedges" << endl;
		vectorPool.clear();
		vector<std::pair<__int64, __int64>> irregulars;
		for (int i = 0; i < val->Vertices.size(); i++)
		{
			vectorPool.push_back(vector<halfedge*>());
		}
		for (auto hf : halfedges)
		{
			vectorPool[hf->P->N].push_back(hf);
		}
		__nVertices = _nVertices;

		for (int i = 0; i < __nVertices; i++)
		{
			auto pool = vectorPool[i];
			vector<__int64> ttt;
			for (auto p : pool)
			{
				ttt.push_back(p->next->P->N);
			}
			std::sort(ttt.begin(), ttt.end());
			vector<__int64> newTTT;
			for (int ii = 0; ii < ttt.size() - 1; ii++)
			{
				if (ttt[ii] == ttt[ii + 1])
				{
					newTTT.push_back(ttt[ii]);
					int tmp = ttt[ii];
					while (ttt[ii] == tmp)ii++;
				}
			}
			for (auto t : newTTT){
				if (t>i)irregulars.push_back(std::make_pair(i, t));
			}
		}
		std::cout <<irregulars.size()<< "irregular halfedges found" << endl;
		for (auto it : irregulars)
		{
			__int64 i;
			__int64 j;
			boost::tie(i, j) = it;
			auto pool1 = vectorPool[i];
			auto pool2 = vectorPool[j];
			vector<halfedge*>go;
			vector<halfedge*>come;
			for (auto s : pool1){
				if (s->next->P->N == j)go.push_back(s);
			}
			for (auto ts : go)
			{
				for (auto s : pool2){
					if (s->pair == ts)come.push_back(s);
				}
			}
			if (go.size() != come.size()){
				std::cout << "error!" << endl;
				std::cin.get();
			}
			if (go.size() >= 2)
			{
				std::cout << "(" << i << "-" << j << ")" << endl;
				__int64 SS = go.size();
				for (__int64 k = 0; k < SS; k++)
				{
					vertices.push_back(new vertex(_nVertices));
					auto newVert = vertices.back();
					auto P = val->Vertices[i];
					auto Q = val->Vertices[j];
					Vector VV(0, 0, 0);
					auto ph = go[k];
					auto cell1 = facet_list[ph->owner->N].first;
					auto D1= facet_list[ph->owner->N].second;
					Point mid((val->Vertices[i].x() + val->Vertices[j].x()) / 2., (val->Vertices[i].y() + val->Vertices[j].y()) / 2., (val->Vertices[i].z() + val->Vertices[j].z()) / 2.);
					auto V1 = cell1->vertex(D1)->point() - mid;
					auto cell2 = facet_list[ph->pair->owner->N].first;
					auto D2 = facet_list[ph->pair->owner->N].second;
					auto V2 = cell2->vertex(D2)->point() - mid;
					VV = VV + V1+V2;
					VV = VV / 2.0/10.;
					vector<halfedge*> pool;
					val->Vertices.push_back(Point((P.x() + Q.x()) / 2. + VV.x(), (P.y() + Q.y()) / 2. + VV.y(), (P.z() + Q.z()) / 2. + VV.z()));
					__int64 l = _nVertices;
					if (i == 6988 && j == 272758)
					{
						std::cout << "(" << i << "-" << l << "-" << j << ")" << endl;
					}
					halfedge* a0, *a1, *b0, *b1;
					auto newHF = new halfedge(newVert);
					auto newHHF = new halfedge(newVert);
					pool.push_back(newHF);
					pool.push_back(newHHF);
					auto newOHHF = new halfedge(ph->prev->P);
					vectorPool[newOHHF->P->N].push_back(newOHHF);
					halfedges.push_back(newHF);
					halfedges.push_back(newHHF);
					halfedges.push_back(newOHHF);
					auto a = ph->prev;
					auto b = ph;
					auto c = newHHF;
					auto A = newOHHF;
					auto B = newHF;
					auto C = ph->next;
					a0 = b;
					a1 = B;
					a->next = b;
					b->next = c;
					c->next = a;
					A->next = B;
					B->next = C;
					C->next = A;
					a->prev = c;
					b->prev = a;
					c->prev = b;
					A->prev = C;
					B->prev = A;
					C->prev = B;
					c->pair = A;
					A->pair = c;
					a->owner = ph->owner;
					b->owner = ph->owner;
					c->owner = ph->owner;
					auto f = new face(_nFaces, ph->owner->corner[0], ph->owner->corner[1], ph->owner->corner[2]);

					for (int ii = 0; ii < 3; ii++)
					{
						if (ph->owner->corner[ii] == j){
							ph->owner->corner[ii] = l;
						}
					}
					val->Faces[ph->owner->N].A = ph->owner->corner[0];
					val->Faces[ph->owner->N].B = ph->owner->corner[1];
					val->Faces[ph->owner->N].C = ph->owner->corner[2];
					if (i == 6988 && j == 272758)
					{
						std::cout << "F:" << val->Faces[ph->owner->N].A << "-" << val->Faces[ph->owner->N].B << "-" << val->Faces[ph->owner->N].C << endl;
					}
					A->owner = f;
					B->owner = f;
					C->owner = f;
					for (int ii = 0; ii < 3; ii++)
					{
						if (f->corner[ii] == i){
							f->corner[ii] = l;
						}
					}
					faces.push_back(f);
					facet_list.push_back(facet_list[ph->owner->N]);
					val->Faces.push_back(MeshFace(f->corner[0], f->corner[1], f->corner[2], val->Faces[ph->owner->N].N));
					if (i == 6988 && j == 272758)
					{
						std::cout << "F:" << val->Faces[f->N].A << "-" << val->Faces[f->N].B << "-" << val->Faces[f->N].C << endl;
					}

					_nFaces++;
					ph = come[k];
					newHF = new halfedge(newVert);
					newHHF = new halfedge(newVert);
					pool.push_back(newHF);
					pool.push_back(newHHF);
					newOHHF = new halfedge(ph->prev->P);
					vectorPool[newOHHF->P->N].push_back(newOHHF);
					halfedges.push_back(newHF);
					halfedges.push_back(newHHF);
					halfedges.push_back(newOHHF);
					a = ph->prev;
					b = ph;
					c = newHHF;
					A = newOHHF;
					B = newHF;
					C = ph->next;
					b0 = b;
					b1 = B;
					a->next = b;
					b->next = c;
					c->next = a;
					A->next = B;
					B->next = C;
					C->next = A;
					a->prev = c;
					b->prev = a;
					c->prev = b;
					A->prev = C;
					B->prev = A;
					C->prev = B;
					c->pair = A;
					A->pair = c;
					a->owner = ph->owner;
					b->owner = ph->owner;
					c->owner = ph->owner;
					f = new face(_nFaces, ph->owner->corner[0], ph->owner->corner[1], ph->owner->corner[2]);
					for (int ii = 0; ii < 3; ii++)
					{
						if (ph->owner->corner[ii] == i){
							ph->owner->corner[ii] = l;
						}
					}
					val->Faces[ph->owner->N].A = ph->owner->corner[0];
					val->Faces[ph->owner->N].B = ph->owner->corner[1];
					val->Faces[ph->owner->N].C = ph->owner->corner[2];
					if (i == 6988 && j == 272758)
					{
						std::cout << "F:" << val->Faces[ph->owner->N].A << "-" << val->Faces[ph->owner->N].B << "-" << val->Faces[ph->owner->N].C << endl;
					}
					A->owner = f;
					B->owner = f;
					C->owner = f;
					for (int ii = 0; ii < 3; ii++)
					{
						if (f->corner[ii] == j){
							f->corner[ii] = l;
						}
					}
					faces.push_back(f);
					facet_list.push_back(facet_list[ph->owner->N]);
					val->Faces.push_back(MeshFace(f->corner[0], f->corner[1], f->corner[2], val->Faces[ph->owner->N].N));
					if (i == 6988 && j == 272758)
					{
						std::cout << "F:" << val->Faces[f->N].A << "-" << val->Faces[f->N].B << "-" << val->Faces[f->N].C << endl;
					}
					_nFaces++;

					vectorPool.push_back(pool);
					_nVertices++;
					a0->pair = b1;
					b1->pair = a0;
					a1->pair = b0;
					b0->pair = a1;
				}
			}

		}
	}
	void  MeshStructure::halfEdgeAdd(vector<face *>f)
	{
		auto _o = orient::unknown;
		vector<face*> faces;
		neighbors.clear();
		for (auto face : f)
		{
			for (int i = 0; i < 3; i++)
			{
				__int64 I = face->corner[i];
				__int64 J = (i == 2) ? face->corner[0] : face->corner[i + 1];
				getfaces(I, J, faces);
				if (faces.size() == 2)
				{
					if (faces[0] == face)
					{
						if (__orientation[faces[1]->N] != orient::unknown)
						{
							_o = __orientation[faces[1]->N] == orient::positive ? orient::negative : orient::positive;
						}
					}
					if (faces[1] == face)
					{
						if (__orientation[faces[0]->N] != orient::unknown)
						{
							_o = __orientation[faces[0]->N] == orient::positive ? orient::negative : orient::positive;
						}
					}
				}
				else
				{
					getfaces(J, I, faces);
					if (faces.size() != 0)
					{
						if (__orientation[faces[0]->N] != orient::unknown)
						{
							_o = __orientation[faces[0]->N];
						}
					}
				}
			}

			__orientation[face->N] = _o == orient::unknown ? orient::positive : _o;
			//std::cout << face->N << endl;

			//register a halfedge
			if (__orientation[face->N] == orient::positive)
			{
				auto he1 = new halfedge(vertices[face->corner[0]]);
				auto he2 = new halfedge(vertices[face->corner[1]]);
				auto he3 = new halfedge(vertices[face->corner[2]]);
				halfedges.push_back(he1);
				halfedges.push_back(he2);
				halfedges.push_back(he3);
				he1->prev = he3; he1->next = he2; he1->owner = face;
				he2->prev = he1; he2->next = he3; he2->owner = face;
				he3->prev = he2; he3->next = he1; he3->owner = face;
				face->firsthalfedge = he1;
			}

			if (__orientation[face->N] == orient::negative)
			{
				auto he1 = new halfedge(vertices[face->corner[2]]);
				auto he2 = new halfedge(vertices[face->corner[1]]);
				auto he3 = new halfedge(vertices[face->corner[0]]);
				halfedges.push_back(he1);
				halfedges.push_back(he2);
				halfedges.push_back(he3);
				he1->prev = he3; he1->next = he2; he1->owner = face;
				he2->prev = he1; he2->next = he3; he2->owner = face;
				he3->prev = he2; he3->next = he1; he3->owner = face;
				face->firsthalfedge = he1;
			}
		}
		for (auto face : f)
		{
			//list up neighbors that are not oriented
			//list up neighbors that are not oriented
			for (int i = 0; i < 3; i++)
			{
				int I = face->corner[i];
				int J = (i == 2) ? face->corner[0] : face->corner[i + 1];
				getfaces(I, J, faces);
				if (faces.size() == 2)
				{
					if (faces[0] == face)
					{
						if (__orientation[faces[1]->N] == orient::unknown)
						{
							neighbors.push_back(faces[1]);// halfEdgeAdd(faces[1]);
						}
					}
					if (faces[1] == face)
					{
						if (__orientation[faces[0]->N] == orient::unknown)
						{
							neighbors.push_back(faces[0]);// halfEdgeAdd(faces[0]);
						}
					}
				}
				else
				{
					getfaces(J, I, faces);
					if (faces.size() != 0)
					{
						if (__orientation[faces[0]->N] == orient::unknown)
						{
							neighbors.push_back(faces[0]);// halfEdgeAdd(faces[0]);
						}
					}
				}
			}
		}
		std::sort(neighbors.begin(),neighbors.end());
		neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());
	}

	void  MeshStructure::halfEdgeAdd(face *f)
	{
		auto _o = orient::unknown;
		vector<face*> faces;
		for (int i = 0; i < 3; i++)
		{
			__int64 I = f->corner[i];
			__int64 J = (i == 2) ? f->corner[0] : f->corner[i + 1];
			getfaces(I, J, faces);
			if (faces.size() == 2)
			{
				if (faces[0] == f)
				{
					if (__orientation[faces[1]->N] != orient::unknown)
					{
						_o = __orientation[faces[1]->N] == orient::positive ? orient::negative : orient::positive;
					}
				}
				if (faces[1] == f)
				{
					if (__orientation[faces[0]->N] != orient::unknown)
					{
						_o = __orientation[faces[0]->N] == orient::positive ? orient::negative : orient::positive;
					}
				}
			}
			else
			{
				getfaces(J, I, faces);
				if (faces.size() != 0)
				{
					if (__orientation[faces[0]->N] != orient::unknown)
					{
						_o = __orientation[faces[0]->N];
					}
				}
			}
		}

		__orientation[f->N] = _o == orient::unknown ? orient::positive : _o;
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
		//list up neighbors that are not oriented
		neighbors.clear();

		for (int i = 0; i < 3; i++)
		{
			int I = f->corner[i];
			int J = (i == 2) ? f->corner[0] : f->corner[i + 1];
			getfaces(I, J, faces);
			if (faces.size() == 2)
			{
				if (faces[0] == f)
				{
					if (__orientation[faces[1]->N] == orient:: unknown)
					{
						neighbors.push_back(faces[1]);// halfEdgeAdd(faces[1]);
					}
				}
				if (faces[1] == f)
				{
					if (__orientation[faces[0]->N] == orient::unknown)
					{
						neighbors.push_back(faces[0]);// halfEdgeAdd(faces[0]);
					}
				}
			}
			else
			{
				getfaces(J, I, faces);
				if (faces.size() != 0)
				{
					if (__orientation[faces[0]->N] == orient::unknown)
					{
						neighbors.push_back(faces[0]);// halfEdgeAdd(faces[0]);
					}
				}
			}
		}
	}
	void MeshStructure::faceTableAdd(int i, int j, face* f)
	{
/*		if (_faceTable.coeff(i, j) == NULL)
		{
			_faceTable.coeffRef(i, j) = new vector<face*>();
		}
		_faceTable.coeffRef(i, j)->push_back(f);*/
		_faceTable[i].push_back(std::make_pair(f, j));

	}
	void MeshStructure::faceTableAdd(face* f)
	{
		for (int i = 0; i < 3; i++)
		{
			__int64 I = f->corner[i];
			__int64 J = 0;
			if (i == 2) J = f->corner[0]; else J = f->corner[i + 1];
			faceTableAdd(I, J, f);
		}
	}
	MeshStructure* MeshStructure::CreateFrom(Mesh *val)
	{
		MeshStructure* ret = new MeshStructure();
		ret->Construct(val);
		return ret;
	}
	Poly MeshStructure::CreateFrom_already_oriented(Mesh *val, vector<Delaunay::Facet> facet_list)
	{
		MeshStructure* ret = new MeshStructure();
		ret->Construct_already_oriented(val, facet_list);
		Poly P;
		Builder<HalfedgeDS> builder;
		builder.mesh = val;
		builder.MS = ret;
		P.delegate(builder);
		delete ret;

		std::cout << "P.is_valid=" << (P.is_valid() ? "true" : "false") << endl;
		std::cout << "P.is_closed=" << (P.is_closed() ? "true" : "false") << endl;
		// This is a stop predicate (defines when the algorithm terminates).
		// In this example, the simplification stops when the number of undirected edges
		// left in the surface mesh drops below the specified number (1000)
		CGAL::Unique_hash_map<edge_descriptor,bool> constraint_hmap(false);
		  Constrained_edge_map constraints_map(constraint_hmap);
		SMS::Constrained_placement<SMS::Midpoint_placement<Poly>,
									Constrained_edge_map > placement(constraints_map);

		// map used to check that constrained_edges and the points of its vertices
		// are preserved at the end of the simplification
		// Warning: the computation of the dihedral angle is only an approximation and can
		//          be far from the real value and could influence the detection of sharp
		//          edges after the simplification
		std::map<edge_descriptor,std::pair<PP, PP> >constrained_edges;
		std::size_t nb_sharp_edges=0;
		edge_iterator eb,ee;
		for(boost::tie(eb,ee) = CGAL::edges(P); eb != ee ; ++eb )
		{
			halfedge_descriptor hd = CGAL::halfedge(*eb,P);
			if ( is_border(*eb,P) ){
				//std::cout << "border" << std::endl;
				++nb_sharp_edges;
				constraint_hmap[*eb]=true;
				constrained_edges[*eb]=std::make_pair(point(source(hd,P),P),
                                            point(target(hd,P),P));
			}else{
				double angle = CGAL::Mesh_3::dihedral_angle(point(target(opposite(hd,P),P),P),
                                                  point(target(hd,P),P),
                                                  point(target(next(hd,P),P),P),
                                                  point(target(next(opposite(hd,P),P),P),P));
				if ( CGAL::abs(angle)<100 ){
					++nb_sharp_edges;
					constraint_hmap[*eb]=true;
					PP p = point(source(hd,P),P);
					PP q = point(target(hd,P),P);
					constrained_edges[*eb]=std::make_pair(p,q);
					//std::cout << "2 " << p << " "  << q << "\n";
				}
			}
		}
		std::cout << "# sharp edges = " << nb_sharp_edges << std::endl;


  		SMS::Count_stop_predicate<Poly> stop(P.size_of_halfedges() / 20);

		// This the actual call to the simplification algorithm.
		// The surface mesh and stop conditions are mandatory arguments.
		// The index maps are needed because the vertices and edges
		// of this surface mesh lack an "id()" field.
		int r = SMS::edge_collapse
			(P
			, stop
			, CGAL::vertex_index_map(get(CGAL::vertex_external_index, P))
			.halfedge_index_map(get(CGAL::halfedge_external_index, P))
			.edge_is_constrained_map(constraints_map)
			.get_placement(SMS::Midpoint_placement<Poly>())
			);

		std::cout << "\nFinished...\n" << r << " edges removed.\n"
			<< (P.size_of_halfedges() / 2) << " final edges.\n";
		std::cout << "P.is_valid=" << (P.is_valid() ? "true" : "false") << endl;
		std::cout << "P.is_closed=" << (P.is_closed() ? "true" : "false") << endl;
		return P;
	}
	void MeshStructure::Clear()
	{
		if (vertices.size() != 0)
		{
			for (auto v : vertices)
			{
				delete(v);
			}
		}
		vertices.clear();
		if (faces.size() != 0)
		{
			for (auto f : faces)
			{
				delete(f);
			}
		}
		faces.clear();
		if (halfedges.size() != 0)
		{
			for (auto hf : halfedges)
			{
				delete(hf);
			}
		}
		/*
		if (_faceTable.nonZeros() != 0)
		{
			for (int k = 0; k < _faceTable.outerSize(); ++k)
			{
				for (SparseMatrix<vector<face*>*>::InnerIterator it(_faceTable, k); it; ++it)
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
		*/
		halfedges.clear();
		innerVertices.clear();
		outerVertices.clear();
		boundaryStart.clear();
	}
	MeshStructure::~MeshStructure()
	{
		this->Clear();
	}
	bool is_border(edge_descriptor e, const Poly& sm)
	{
		return (CGAL::face(CGAL::halfedge(e, sm), sm) == boost::graph_traits<Poly>::null_face())
			|| (CGAL::face(opposite(CGAL::halfedge(e, sm), sm), sm) == boost::graph_traits<Poly>::null_face());
	}
	PP point(vertex_descriptor vd, const Poly& sm)
	{
		return get(CGAL::vertex_point, sm, vd);
	}
}
