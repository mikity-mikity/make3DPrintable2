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
        int N;
        int corner[3];
        halfedge *firsthalfedge;
		face();
		void face::setup(int _N, int A, int B, int C);
		face(int _N, int A, int B, int C);
    };
}
#endif