#include"vertex.h"
#include"halfedge.h"
namespace GeometryProcessing
{
	vertex::vertex(int _N)
	{
		N = _N;
		hf_begin=NULL;
	}
	bool vertex::isInner()
	{
		return onering[0] == onering[onering.size() - 1]->next;
	}
	bool vertex::isBoundary()
	{
		return onering[0] != onering[onering.size() - 1]->next;
	}
}