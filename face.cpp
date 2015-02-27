#include"face.h"
namespace GeometryProcessing{
	face::face()
	{
		firsthalfedge = NULL;
	}
	void face::setup(int _N, int A, int B, int C)
	{
		corner[0] = A;
		corner[1] = B;
		corner[2] = C;
		N = _N;
	}
	face::face(int _N, int A,int B,int C)
	{
		corner[0]=A;
		corner[1]=B;
		corner[2]=C;
		N = _N;
		firsthalfedge=NULL;
	}
}
