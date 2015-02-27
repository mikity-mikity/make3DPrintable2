#include "halfedge.h"
#include "vertex.h"
namespace GeometryProcessing{

    //Warning
    //After boundary halfedges are appropriately setup, Naked halfedges should not exist.
	bool halfedge::isNaked()
	{
		if(pair==NULL)return true;
		return false;
	}
	bool halfedge::isBoundary()
	{
		if(owner==NULL)return true;
		return false;
	}
	bool halfedge::ifPairIsBoundary()
	{
		if(pair->owner==NULL)return true;
		return false;
	}
	void halfedge::setVertex(vertex *_P)
	{
		P = _P;
		if (_P->hf_begin == NULL) _P->hf_begin = this;
	}
	halfedge::halfedge()
	{
		P = NULL;
		owner = NULL;
		pair = NULL;
		next = NULL;
		prev = NULL;
	}
	halfedge::halfedge(vertex *_P)
	{
		P = _P;
		if (_P->hf_begin == NULL) _P->hf_begin = this;
		owner=NULL;
		pair=NULL;
		next=NULL;
		prev=NULL;
	}
}
