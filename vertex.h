#ifndef _VERTEX_H_
#define _VERTEX_H_
#include <vector>
#include "halfedge.h"
using namespace std;

namespace GeometryProcessing
{
	class halfedge;
	class vertex
    {
	public:
		int N;
        vector<halfedge*> star;
        vector<halfedge*> onering;
        halfedge* hf_begin;
        vertex(int _N);
        bool isInner();
        bool isBoundary();
    };
}
#endif
