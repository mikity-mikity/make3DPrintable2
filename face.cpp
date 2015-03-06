#include"face.h"
namespace GeometryProcessing{
	face::face()
	{
		firsthalfedge = NULL;
	}
	void face::setup(__int64 _N, __int64 A, __int64 B, __int64 C)
	{
		corner[0] = A;
		corner[1] = B;
		corner[2] = C;
		N = _N;
	}
	face::face(__int64 _N, __int64 A,__int64 B,__int64 C)
	{
		corner[0]=A;
		corner[1]=B;
		corner[2]=C;
		N = _N;
		firsthalfedge=NULL;
	}
}
