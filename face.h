#ifndef _FACE_H_
#define _FACE_H_
#include <vector>
#include "halfedge.h"
using namespace std;
namespace GeometryProcessing
{
	class halfedge;
	class face
    {
	public:
        __int64 N;
        __int64 corner[3];
        halfedge *firsthalfedge;
		face();
		void face::setup(__int64 _N, __int64 A, __int64 B, __int64 C);
		face(__int64 _N, __int64 A, __int64 B, __int64 C);
    };
}
#endif