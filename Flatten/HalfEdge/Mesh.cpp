/** \file Mesh.cpp Definitions for the Mesh class.
    \author Liliya Kharevych
    \date March 2006
*/

#include "Mesh.h"
#include <iostream>
#include <fstream>
#include "util.h"
#include <stack>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>
#include <assert.h>

using namespace std;

Mesh::Mesh(void)
{
   numBoundaryVertices = 0;
}

Mesh::~Mesh(void)
{
   clear();
}

void Mesh::clear() {
   VertexIterator vIter = vertexIterator();
   FaceIterator fIter = faceIterator();
   HalfEdgeIterator eIter = halfEdgeIterator();

   while (!vIter.end()) {
      delete vIter.vertex();
      vIter++;
   }
   while (!fIter.end()) {
      delete fIter.face();
      fIter++;
   }
   while (!eIter.end()) {
      delete eIter.half_edge();
      eIter++;
   }
}

Vertex * Mesh::addVertex(const Eigen::Vector3d & _p) {
   Vertex * v = new Vertex();
   v->p = _p;
   v->ID = (int)vertices.size();
   vertices.push_back(v);
   return v;
}

Face * Mesh::addFace(std::vector<int> faceVerts) {
   Face * f = new Face();
   f->ID = (int)faces.size();

   Edge * firstEdge = NULL;
   Edge * last = NULL;
   Edge * current = NULL;

   unsigned int i;

   for (i = 0; i < faceVerts.size()-1; i++) {
      current = addEdge(faceVerts[i], faceVerts[i+1]);
   
      check_error (!current, "Problem while parsing mesh faces.");

      current->face = f;

      if (last)
	 last->next = current;
      else
	 firstEdge = current;
      last = current;
   }

   current = addEdge (faceVerts[i], faceVerts[0]);
   check_error (!current, "Problem while parsing mesh faces.");
   
   current->face = f;

   last->next = current;
   current->next = firstEdge;

   f->edge = firstEdge;
   faces.push_back(f);
      
   return f;
}

Edge * Mesh::addEdge (int i, int j) {
   Edge eTmp;
   eTmp.vertex = vertices[i];

   Edge eTmpPair;
   eTmpPair.vertex = vertices[j];
   eTmp.pair = &eTmpPair;

   Mesh::EdgeSet::iterator eIn = edges.find(&eTmp);
   Edge * e;

   if (eIn != edges.end()){
      e = *eIn;
      if (e->face != NULL)
	 return NULL;
   }
   else {
      e = new Edge();
      Edge * pair = new Edge();

      e->vertex = vertices[i];
      pair->vertex = vertices[j];

      pair->pair = e;
      e->pair = pair;

      edges.insert(e);
      edges.insert(pair);

      pair->vertex->edge = pair;
   }   
   return e;
}

void Mesh::computeNormals()
{
	map< int, Eigen::Vector3d >					normals;

	for( Mesh::FaceIterator fIter = this->faceIterator() ; !fIter.end() ; fIter++ )
	{
		Face		*face = fIter.face();
		Eigen::Vector3d		u = face->edge->next->vertex->p - face->edge->vertex->p;
		Eigen::Vector3d		v = face->edge->next->next->vertex->p - face->edge->vertex->p;
		Eigen::Vector3d		n = u.cross(v).normalized();
		normals.insert( pair<int, Eigen::Vector3d>(face->ID, n) );
	}

	for( Mesh::VertexIterator vIter = this->vertexIterator() ; !vIter.end() ; vIter++ )
	{
		int		valence = 0;
		vIter.vertex()->n = Eigen::Vector3d(0,0,0);

		for( Vertex::EdgeAroundIterator	around = vIter.vertex()->iterator(); !around.end(); around++ )
		{
			if(around.edge_out()->face)
			{
				int		faceId = around.edge_out()->face->ID;
				vIter.vertex()->n += normals[faceId];
				valence++;
			}
		}

		vIter.vertex()->n /= (double)valence;
		vIter.vertex()->n.normalize();
	}
}


/** Called after the mesh is created to link boundary edges */
void Mesh::linkBoundary() {
   HalfEdgeIterator hIter = halfEdgeIterator();

   for (; !hIter.end(); hIter++) {
      if (!hIter.half_edge()->face) 
	 hIter.half_edge()->vertex->edge = hIter.half_edge();
   }

    VertexIterator vIter = vertexIterator();
   Edge *next, *beg;
   while (!vIter.end()) {
      if (vIter.vertex()->isBoundary()) {
         beg = vIter.vertex()->edge;
         next = beg;
         while (beg->pair->face)
            beg = beg->pair->next;
         beg = beg->pair;
         beg->next = next;
         numBoundaryVertices++;
      }
      vIter++;
   }

}

/** Computes information about the mesh:

Number of boundary loops,
Number of connected components,
Genus of the mesh.

Only one connected compoment with 0 or 1 boundary loops can be parameterize.
Singularities should be assigned in such a way that Gauss-Bonet theorem is satisfied.
*/
void Mesh::computeMeshInfo() {
   cout << "Topological information about the mesh:" << endl;
   // Number of boundary loops
   Mesh::HalfEdgeIterator hIter = halfEdgeIterator();
   for (; !hIter.end(); hIter++) {
      hIter.half_edge()->check = false;
   }
   numBoundaryLoops = 0;
   for (hIter.reset(); !hIter.end(); hIter++) {
      Edge * e = hIter.half_edge();
      if (e->face)
	 e->check = true;
      else if (!e->check) {
	 Edge * beg = NULL;
	 while (e != beg) {
	    if (!beg) beg = e;
	    check_error(!e, "Building the mesh failed, probem with input format.");
	    e->check = true;
	    e = e->next;
	 }
	 numBoundaryLoops++;
      }
   }
   cout << "Mesh has " << numBoundaryLoops << " boundary loops." << endl;
   // Number of connected components
   numConnectedComponents = 0;
   Mesh::FaceIterator fIter = faceIterator();
   for (; !fIter.end(); fIter++) {
      fIter.face()->check = false;
   }
   stack<Edge *> toVisit;
   for (fIter.reset(); !fIter.end(); fIter++) {
      if (!fIter.face()->check) {
	 numConnectedComponents++;
	 toVisit.push(fIter.face()->edge);
	 while (!toVisit.empty()) {
	    Face * fIn = toVisit.top()->face; 
	    toVisit.pop();
	    fIn->check = true;     
	    Face::EdgeAroundIterator iter = fIn->iterator();
	    for (; !iter.end(); iter++) 
	       if (iter.edge()->pair->face && !iter.edge()->pair->face->check)
		  toVisit.push(iter.edge()->pair);
	 }
      }
   }
   cout << "Mesh has " << numConnectedComponents << " connected components." << endl;
   // Genus number
   check_error(numConnectedComponents == 0, "The mesh read is empty.");
   numGenus = 
      (1 - (numVertices() - numEdges() + numFaces() + numBoundaryLoops ) / 2) / numConnectedComponents;
   cout << "Mesh is genus " << numGenus << "." << endl;
}

/** Check if all the vertices in the mesh have at least on edge coming out of them */
bool Mesh::checkVertexConection() {
   Mesh::FaceIterator fIter = faceIterator();
   Mesh::VertexIterator vIter = vertexIterator();
   bool conectedVert = true;

   for (;!vIter.end(); vIter++)
      vIter.vertex()->check = false;

   for (fIter.reset(); !fIter.end(); fIter++) {
      Face::EdgeAroundIterator around = fIter.face()->iterator();
      for (;!around.end();around++)
	 around.vertex()->check = true;
   }
   for (vIter.reset(); !vIter.end(); vIter++) {
      if (!vIter.vertex()->check) {
	 cerr << "Vertex " << vIter.vertex()->ID << " is not connected." << endl;
	 conectedVert = false;
      }
   }

   return conectedVert;
}

/** Manifoldness check: only one disk should be adjusted on any vertex */
bool Mesh::checkManifold() {
   Mesh::HalfEdgeIterator eIter = halfEdgeIterator();
   Mesh::VertexIterator vIter = vertexIterator();
   bool manifold = true;

   for (;!eIter.end(); eIter++)
      eIter.half_edge()->check = false;

   for (vIter.reset(); !vIter.end(); vIter++) {
      Vertex::EdgeAroundIterator around = vIter.vertex()->iterator();
      for (;!around.end();around++)
	 around.edge_out()->check = true;
   }

   for (eIter.reset(); !eIter.end(); eIter++) {
      if (!eIter.half_edge()->check) {
	 cerr << "Mesh is not manifold - more then one disk at vertex " 
	    << eIter.half_edge()->vertex->ID << endl;
	 manifold = false;
	 break;
      }
   }

   return manifold;
}

void Mesh::CenterAndNormalize()
{
	double maxX, maxY, maxZ, minX, minY, minZ;
	maxX = maxY = maxZ = -1e38;
	minX = minY = minZ = 1e38;

	for( VertexIterator vIter = vertexIterator() ; !vIter.end() ; vIter++ )
	{
		if(vIter.vertex()->p.x() > maxX) maxX = vIter.vertex()->p.x();
		if(vIter.vertex()->p.y() > maxY) maxY = vIter.vertex()->p.y();
		if(vIter.vertex()->p.z() > maxZ) maxZ = vIter.vertex()->p.z();
		if(vIter.vertex()->p.x() < minX) minX = vIter.vertex()->p.x();
		if(vIter.vertex()->p.y() < minY) minY = vIter.vertex()->p.y();
		if(vIter.vertex()->p.z() < minZ) minZ = vIter.vertex()->p.z();
	}

	Eigen::Vector3d min(minX,minY,minZ);
	Eigen::Vector3d max(maxX,maxY,maxZ);

	Eigen::Vector3d center = min + (max - min) * 0.5;

	for( VertexIterator vIter = vertexIterator() ; !vIter.end() ; vIter++ )
	{
		vIter.vertex()->p -= center;
	}

	double diag = (max - min).norm() / 2.0;
	double scale = 1.0 / diag;
	for( VertexIterator vIter = vertexIterator() ; !vIter.end() ; vIter++ )
	{
		vIter.vertex()->p *= scale;
	}
}


/** Loads mesh from obj file

Standard format for OBJ file is

v double double double

v double double double

f int int int

f int int int

Files with vt tags also can be parsed
*/
void Mesh::readOBJ(const char * obj_file) {
   string front;
   string v = "v", vt = "vt", f = "f";
   Eigen::Vector3d vert;
   vector<int> verts;
   vector<Eigen::Vector3d> uvVec;
   vector<int> uvs;
   char etc;
   int id;

   ifstream in(obj_file);

   check_error(!in, "Cannot find input obj file.");

   bool hasUV = false;

   while(!in.eof() || in.peek() != EOF) {
      in >> ws; //extract whitespaces
      if (in.eof() || in.peek() == EOF)
	 break;
      if (in.peek() == '#') {
	 in.ignore(300, '\n');
      }
      else {
	 in >> front;
	 if (front == v) {
	    in >> vert.x() >> vert.y() >> vert.z();
	    addVertex(vert);
	 }
	 else if (front == vt){
	    in >> vert.x() >> vert.y();
	    vert.z() = 0;
	    uvVec.push_back(vert);
	    hasUV = true;
	 }
	 else if (front == f) {
	    verts.clear();
	    uvs.clear();
	    while (in >> id) {

	       check_error(id > numVertices(), "Problem with input OBJ file.");

	       verts.push_back(id-1);
	       bool vtNow = true;
	       if (in.peek() == '/'){
    		  in >> etc;
		  if (in.peek() != '/') {
      		     in >> id;
		     check_warn(id > numVertices(), "Texture coordinate index is greater then number of vertices.");
		     if (id < numVertices() && hasUV) {
			uvs.push_back(id-1);
			vtNow = false;
		     }
		  }
	       }
	       if (in.peek() == '/'){
		  int tmp;
   		  in >> etc;
		  in >> tmp;
	       }
	       if (hasUV && vtNow) {
		  uvs.push_back(id-1);
	       }
	    }
	    in.clear(in.rdstate() & ~ios::failbit);
	    Face * f = addFace(verts);

	    if (hasUV && uvs.size() != 0){
	       int k = 0;
	       for (Face::EdgeAroundIterator e = f->iterator(); !e.end(); e++, k++)
		  e.vertex()->uv = uvVec[uvs[k]];
	    }
	 }
	 else {
	    string line;
	    getline(in, line);
	    cout << "Warning, line: " << line << " ignored." << endl;	
	 }
      }
   }

   in.close();

   // Finnish building the mesh, should be called after each parse.
   finishMesh();
}

/* Write mesh in OBJ format to obj_file parameter */
void Mesh::writeOBJ(const char * obj_file) {
   ofstream out(obj_file);
   check_error (!out, "Cannot find output file.");

   Mesh::VertexIterator vItr = vertexIterator();
   for (vItr.reset(); !vItr.end(); vItr++)
      out << "v " << vItr.vertex()->p << endl;

   for (vItr.reset(); !vItr.end(); vItr++)
      out << "vt " << vItr.vertex()->uv.x() << " " << vItr.vertex()->uv.y() << endl;

   Mesh::FaceIterator fItr = faceIterator();
   for (fItr.reset(); !fItr.end(); fItr++) {
      Face::EdgeAroundIterator around = fItr.face()->iterator();
      out << "f";
      for ( ; !around.end(); around++)
	 out << " " << (around.vertex()->ID + 1) << "/" << (around.vertex()->ID + 1);
      out << endl;
   }
}

/* Write mesh in VT format (only text coodrinates) to vt_file parameter */
void Mesh::writeVT(const char * vt_file) {
   ofstream out(vt_file);
   check_error (!out, "Cannot find output file.");

   Mesh::VertexIterator vItr = vertexIterator();
   for (vItr.reset(); !vItr.end(); vItr++)
      out << "vt " << vItr.vertex()->uv.x() << " " << vItr.vertex()->uv.y() << endl;
}
