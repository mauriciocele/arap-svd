// Flatten.cpp : Defines the entry point for the console application.
//
#ifdef WIN32
#define NOMINMAX
#include <windows.h>
#endif

#if defined (__APPLE__) || defined (OSX)
	#include <OpenGL/gl.h>
	#include <GLUT/glut.h>
#else
	#include <GL/gl.h>
	#include <GL/glut.h>
#endif

#include "GA/c3ga.h"
#include "GA/c3ga_util.h"
#include "GA/gl_util.h"

#include "primitivedraw.h"
#include "gahelper.h"
#include "Laplacian.h"

#include <memory>

#include <vector>
#include <map>
#include "numerics.h"
#include "HalfEdge/Mesh.h"
#include "GARotorEstimator.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

const char *WINDOW_TITLE = "Interactive 3D Shape Deformation using Conformal Geometric Algebra";

// GLUT state information
int g_viewportWidth = 800;
int g_viewportHeight = 600;

void display();
void reshape(GLint width, GLint height);
void MouseButton(int button, int state, int x, int y);
void MouseMotion(int x, int y);
void KeyboardUpFunc(unsigned char key, int x, int y);
void SpecialFunc(int key, int x, int y);
void SpecialUpFunc(int key, int x, int y);
void Idle();
void DestroyWindow();

using namespace c3ga;
using namespace std;
using namespace numerics;

class Camera
{
public:
	float		pos[3];
	float		fw[3];
	float		up[3];
	float		translateVel;
	float		rotateVel;

	Camera()
	{
		float		_pos[] = { 0, 0, -2};
		float		_fw[] = { 0, 0, 1 };
		float		_up[] = { 0, 1, 0 };

		translateVel = 0.005;
		rotateVel = 0.005;
		memcpy(pos, _pos, sizeof(float)*3);
		memcpy(fw, _fw, sizeof(float)*3);
		memcpy(up, _up, sizeof(float)*3);
	}

	void glLookAt()
	{
		gluLookAt( pos[0], pos[1], pos[2], fw[0],  fw[1],  fw[2], up[0],  up[1],  up[2] );
	}
};

class VertexBuffer
{
public:
	std::vector<Eigen::Vector3d> deformedPosition; //deformed mesh position
	std::vector<normalizedPoint> position; //primary backup
	std::vector<normalizedPoint> positionConstrained; //primary backup
	std::vector<Eigen::Vector3d> laplacianCoordinate; //laplacian Coordinate
	std::vector<Eigen::Vector3d> normal; //for rendering (lighting)
	std::vector<Eigen::Vector3d> normalOrig; //primary backup
	std::vector<Eigen::Quaterniond> R;
	std::vector<bool> isBoundary;
	int size;

	VertexBuffer() : size(0)
	{
	}

	void resize(int size)
	{
		this->size = size;
		deformedPosition.resize(size);
		position.resize(size);
		positionConstrained.resize(size);
		laplacianCoordinate.resize(size);
		normal.resize(size);
		normalOrig.resize(size);
		R.resize(size);
		isBoundary.resize(size);
	}
	int get_size() { return size; }
};

class IndexBuffer{
public:
	std::vector<int> faces;
	int size;

	IndexBuffer() : size(0)
	{
	}

	void resize(int size)
	{
		this->size = size;
		faces.resize(size);
	}
	int get_size() { return size; }

};

class Handle
{
public:
	TRversor R;
	TRversor T;
	dualSphere dS;
	bool fixed;
	dualPlane Plane;

	Handle(dualSphere dS, bool fixed, dualPlane dP)
	{
		normalizedPoint x = DualSphereCenter(dS);
		T = _TRversor( exp(-0.5 * _vectorE3GA(x)^ni) );
		R = _TRversor(1.0);
		Plane = dP;
		this->dS = inverse(T) * dS * T;
		this->fixed = fixed;
	}

	TRversor GetTRVersor()
	{
		return T * R;
	}
};

Camera g_camera;
Mesh mesh;
vectorE3GA g_prevMousePos;
bool g_rotateModel = false;
bool g_rotateModelOutOfPlane = false;
rotor g_modelRotor = _rotor(1.0);
bool g_rotateKeyRotors = false;
bool g_translateKeyRotors = false;
bool g_computeBasis = false;
float g_dragDistance = -1.0f;
int g_dragObject;
std::shared_ptr<SparseMatrix> A;
std::set<int> constraints;
std::set<int> fconstraints;
Eigen::SparseMatrix<double> Lc;
Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
int systemType = LaplaceBeltrami; //MeanValue; //LaplaceBeltrami
bool g_showSpheres = true;
bool g_showWires = false;
bool g_iterateManyTimes = false;
Eigen::MatrixXd b3;
Eigen::MatrixXd xyz;

VertexBuffer vertexDescriptors;
IndexBuffer triangles;
std::vector<std::shared_ptr<Handle>> handles;

Eigen::Vector3d _vector3d( const normalizedPoint& p)
{
	return Eigen::Vector3d(p.e1(), p.e2(), p.e3());
}

std::vector<Eigen::Vector3d> ComputeLaplacianCoordinates(std::shared_ptr<SparseMatrix> A, Mesh* mesh)
{
	std::vector<Eigen::Vector3d> laplacianCoordinates(A->numColumns(), Eigen::Vector3d(0,0,0));
	
	auto numRows = A->numRows();

	for( int i = 0; i < numRows ; ++i)
	{
		SparseMatrix::RowIterator aIter = A->iterator(i);
		for( ; !aIter.end() ; ++aIter )
		{
			auto j = aIter.columnIndex();
			laplacianCoordinates[i] += mesh->vertexAt(j)->p * aIter.value();
		}
	}
	
	return laplacianCoordinates;
}

void PreFactor(std::shared_ptr<SparseMatrix> A)
{
	Lc = Eigen::SparseMatrix<double>(A->numRows(), A->numColumns());

	auto numRows = A->numRows();
	for( int i = 0; i < numRows ; ++i)
	{
		if(constraints.find(i) == constraints.end() && fconstraints.find(i) == fconstraints.end())
		{
			SparseMatrix::RowIterator aIter = A->iterator(i);
			for( ; !aIter.end() ; ++aIter )
			{
				auto j = aIter.columnIndex();
				Lc.insert(i, j) = (*A)(i,j);
			}
		}
		else
		{
			Lc.insert(i, i) = 1.0;
		}
	}
	Lc.makeCompressed();
	solver.compute(Lc);
	if(solver.info() != Eigen::Success) {
		// TODO: error handling
	}
}

int main(int argc, char* argv[])
{
	mesh.readOBJ("cactus2.obj"); //armadillo-5k-smooth.obj female.obj david.obj rabbit.obj tyra.obj horse.obj cylinder.obj bar.obj planewithpeaks.obj dragon.obj catHead.obj  cactus.obj  bunny.obj  764_hand-olivier-10kf.obj armadillo.obj

	mesh.CenterAndNormalize();

	mesh.computeNormals();

	// GLUT Window Initialization:
	glutInit (&argc, argv);
	glutInitWindowSize(g_viewportWidth, g_viewportHeight);
	glutInitDisplayMode( GLUT_RGB | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow(WINDOW_TITLE);

	// Register callbacks:
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMouseFunc(MouseButton);
	glutMotionFunc(MouseMotion);
	glutKeyboardUpFunc(KeyboardUpFunc);
	glutSpecialFunc(SpecialFunc);
	glutSpecialUpFunc(SpecialUpFunc);
	glutIdleFunc(Idle);
	atexit(DestroyWindow);

	InitializeDrawing();

	dualPlane P1 = _dualPlane(c3gaPoint(.0, 0.8,.0) << (_vectorE3GA(0,1,0)*ni));
	dualPlane P2 = _dualPlane(c3gaPoint(.0,-0.8,.0) << (_vectorE3GA(0,-1,0)*ni));

	//cactus.obj
	handles.push_back(std::shared_ptr<Handle>(new Handle(_dualSphere(c3gaPoint(.0, .04, 0.7) - 0.5*SQR(0.15)*ni), false, P1)));
	handles.push_back(std::shared_ptr<Handle>(new Handle(_dualSphere(c3gaPoint(.0, .04, -0.8) - 0.5*SQR(0.15)*ni), true, P2)));

	////cylinder.obj
	//handles.push_back( boost::shared_ptr<Handle>( new Handle( _dualSphere(c3gaPoint(.0, 0.9,.0) - 0.5*SQR(0.25)*ni), false, P1 ) ) );
	//handles.push_back( boost::shared_ptr<Handle>( new Handle( _dualSphere(c3gaPoint(.0,-0.9,.0) - 0.5*SQR(0.25)*ni), true, P2 ) ) );

	//Armadillo Pie y Mano
	//handles.push_back( boost::shared_ptr<Handle>( new Handle( _dualSphere(c3gaPoint(-.5, 0.45,-.3) - 0.5*SQR(0.15)*ni), false, P1 ) ) );
	//handles.push_back( boost::shared_ptr<Handle>( new Handle( _dualSphere(c3gaPoint(.2,-0.6,.1) - 0.5*SQR(0.15)*ni), true, P2 ) ) );
	//Armadillo Pubis y Cabeza
	//handles.push_back( boost::shared_ptr<Handle>( new Handle( _dualSphere(c3gaPoint(.0, 0.4,-.2) - 0.5*SQR(0.15)*ni), false, P1 ) ) );
	//handles.push_back( boost::shared_ptr<Handle>( new Handle( _dualSphere(c3gaPoint(.0,-0.05,.1) - 0.5*SQR(0.15)*ni), true, P2 ) ) );

	//handles[0]->extrema = 0;
	//handles[1]->extrema = 1;
	//extremas[0]->handle = 0;
	//extremas[1]->handle = 1;

	vertexDescriptors.resize(mesh.numVertices());
	for( Mesh::VertexIterator vIter = mesh.vertexIterator() ; !vIter.end() ; vIter++ )
	{
		int i = vIter.vertex()->ID;
		vertexDescriptors.position[i] = c3gaPoint( vIter.vertex()->p.x(), vIter.vertex()->p.y(), vIter.vertex()->p.z() );
		vertexDescriptors.normalOrig[i] = vIter.vertex()->n;
		vertexDescriptors.isBoundary[i] = vIter.vertex()->isBoundary();

		for( int h = 0 ; h < handles.size() ;  ++h )
		{
			TRversor TR = handles[h]->GetTRVersor();

			if( _double(vertexDescriptors.position[i] << (TR * handles[h]->dS * inverse(TR))) > 0 ) //inside the sphere
			//if( _double(vd->position << handles[i]->Plane) > 0 ) //positive side of the plane
			{
				if( !handles[h]->fixed )
				{
					constraints.insert(i);
				}
				else
				{
					fconstraints.insert(i);
				}
			}
		}
		if(vertexDescriptors.isBoundary[i])
		{
			fconstraints.insert(i);
		}
	}

	triangles.resize(mesh.numFaces()*3);
	for( Mesh::FaceIterator fIter = mesh.faceIterator() ; !fIter.end() ; fIter++ )
	{
		Face		*face = fIter.face();
		int i = face->ID;
		int	v1 = face->edge->vertex->ID;
		int	v2 = face->edge->next->vertex->ID;
		int	v3 = face->edge->next->next->vertex->ID;
		triangles.faces[i*3 + 0] = v1;
		triangles.faces[i*3 + 1] = v2;
		triangles.faces[i*3 + 2] = v3;
	}

	A = CreateLaplacianMatrix( &mesh, systemType );
	
	vertexDescriptors.laplacianCoordinate = ComputeLaplacianCoordinates(A, &mesh);

	TRversor R1 = handles[0]->GetTRVersor();
	for(std::set<int>::iterator citer = constraints.begin(); citer != constraints.end() ; citer++)
		vertexDescriptors.positionConstrained[*citer] = normalize(_point(inverse(R1) * vertexDescriptors.position[*citer] * R1));

	TRversor R2 = handles[1]->GetTRVersor();
	for(std::set<int>::iterator citer = fconstraints.begin(); citer != fconstraints.end() ; citer++)
		vertexDescriptors.positionConstrained[*citer] = normalize(_point(inverse(R2) * vertexDescriptors.position[*citer] * R2));

	b3 = Eigen::MatrixXd(A->numRows(), 3);

	PreFactor(A);

	glutMainLoop();

	return 0;
}

void SolveLinearSystem(VertexBuffer& vertexDescriptors)
{
	int n = vertexDescriptors.get_size();

	for( int i = 0 ; i < n ; ++i )
	{
		b3.row(i) = vertexDescriptors.laplacianCoordinate[i];
	}

	TRversor R2 = handles[1]->GetTRVersor();
	TRversor R2inv = inverse(R2);
	for(std::set<int>::iterator citer = fconstraints.begin(); citer != fconstraints.end() ; citer++)
	{
		int i = *citer;
		auto constraint = normalize(_point(R2 * vertexDescriptors.positionConstrained[i] * R2inv));
		b3(i,0) = constraint.e1();
		b3(i,1) = constraint.e2();
		b3(i,2) = constraint.e3();
	}
	TRversor R1 = handles[0]->GetTRVersor();
	TRversor R1inv = inverse(R1);
	for(std::set<int>::iterator citer = constraints.begin(); citer != constraints.end() ; citer++)
	{
		int i = *citer;
		auto constraint = normalize(_point(R1 * vertexDescriptors.positionConstrained[i] * R1inv));
		b3(i,0) = constraint.e1();
		b3(i,1) = constraint.e2();
		b3(i,2) = constraint.e3();
	}

	xyz = solver.solve(b3);

	for( int i = 0 ; i < n ; ++i )
	{
		vertexDescriptors.deformedPosition[i] = xyz.row(i);
		vertexDescriptors.normal[i] = vertexDescriptors.normalOrig[i];
	}
}

void UpdateLaplaciansRotation(Mesh *mesh, std::shared_ptr<SparseMatrix> A, VertexBuffer& vertexDescriptors)
{
	Eigen::Matrix3d m;
	for( Mesh::VertexIterator vIter = mesh->vertexIterator() ; !vIter.end() ; vIter++ )
	{
		m.setZero();
		int i = vIter.vertex()->ID;
		const Eigen::Vector3d &pi = mesh->vertexAt(i)->p;
		const Eigen::Vector3d &tpi = vertexDescriptors.deformedPosition[i];
		double S = 0;
		for(Vertex::EdgeAroundIterator edgeAroundIter = vIter.vertex()->iterator() ; !edgeAroundIter.end() ; edgeAroundIter++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;
			const Eigen::Vector3d &pj = mesh->vertexAt(j)->p;
			const Eigen::Vector3d &tpj = vertexDescriptors.deformedPosition[j];
			Eigen::Vector3d eij = (pj - pi) * (*A)(i, j);
			Eigen::Vector3d teij = tpj - tpi;
			m += eij * teij.transpose();
			S += eij.dot(eij);
			S += teij.dot(teij) * (*A)(i, j);
		}
		vertexDescriptors.R[i] = GARotorEstimator(m, S);
	}

	for( Mesh::VertexIterator vIter = mesh->vertexIterator() ; !vIter.end() ; vIter++ )
	{
		int i = vIter.vertex()->ID;
		vertexDescriptors.laplacianCoordinate[i].setZero();
		for(Vertex::EdgeAroundIterator edgeAroundIter = vIter.vertex()->iterator() ; !edgeAroundIter.end() ; edgeAroundIter++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;
			double wij = (*A)(i, j);
			Eigen::Quaterniond &Ri = vertexDescriptors.R[i];
			Eigen::Quaterniond &Rj = vertexDescriptors.R[j];
			Eigen::Vector3d V = mesh->vertexAt(j)->p - mesh->vertexAt(i)->p;
			vertexDescriptors.laplacianCoordinate[i] += 0.5 * wij * (Ri._transformVector(V) + Rj._transformVector(V));
		}
	}
}

vectorE3GA Color( double d )
{
	static vectorE3GA	c0 = _vectorE3GA( 1, 1, 1);
	static vectorE3GA	c1 = _vectorE3GA( 1, 1, 0);
	static vectorE3GA	c2 = _vectorE3GA( 0, 1, 0);
	static vectorE3GA	c3 = _vectorE3GA( 0, 1, 1);
	static vectorE3GA	c4 = _vectorE3GA( 0, 0, 1);

	if( d < 0.25 )
	{
		double alpha = (d - 0.0) / (0.25-0.0);
		return (1.0 - alpha) * c0 + alpha * c1;
	}
	else if( d < 0.5 )
	{
		double alpha = (d - 0.25) / (0.5-0.25);
		return (1.0 - alpha) * c1 + alpha * c2;
	}
	else if( d < 0.75 )
	{
		double alpha = (d - 0.5) / (0.75-0.5);
		return (1.0 - alpha) * c2 + alpha * c3;
	}
	else
	{
		double alpha = (d - 0.75) / (1.0-0.75);
		return (1.0 - alpha) * c3 + alpha * c4;
	}
}

void display()
{
	/*
	 *	matrices
	 */
	glViewport( 0, 0, g_viewportWidth, g_viewportHeight );
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	pickLoadMatrix();
	GLpick::g_frustumFar = 1000.0;
	GLpick::g_frustumNear = .1;
	gluPerspective( 60.0, (double)g_viewportWidth/(double)g_viewportHeight, GLpick::g_frustumNear, GLpick::g_frustumFar );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glShadeModel(GL_SMOOTH);	//gouraud shading
	glClearDepth(1.0f);
	glClearColor( .75f, .75f, .75f, .0f );
	glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );

	/*
	 *	estados
	 */
	glEnable(GL_CULL_FACE);		//face culling
	glCullFace( GL_BACK );
	glFrontFace( GL_CCW );
	glEnable(GL_DEPTH_TEST);	//z-buffer
	glDepthFunc(GL_LEQUAL);

	/*
	 *	iluminacion
	 */
	float		ambient[] = { .3f, .3f, .3f, 1.f };
	float		diffuse[] = { .3f, .3f, .3f, 1.f };
	float		position[] = { .0f, 0.f, -150.f, 1.f };
	float		specular[] = { 1.f, 1.f, 1.f };

	glLightfv( GL_LIGHT0, GL_AMBIENT, ambient );
	glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuse );
	glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0);
	glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.0125);
	glEnable(  GL_LIGHT0   );
	glEnable(  GL_LIGHTING );
	glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specular );
	glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, 50.f );

	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glLoadIdentity();

	g_camera.glLookAt();

	glLightfv( GL_LIGHT0, /*GL_POSITION*/GL_SPOT_DIRECTION, position );
	glPushMatrix();

	rotorGLMult(g_modelRotor);

	if(!g_computeBasis)

	{
		static bool oneTime = true;
		if(oneTime || (g_rotateKeyRotors || g_translateKeyRotors))
		{
			if(oneTime == true)
			{
				SolveLinearSystem(vertexDescriptors);
			}

			oneTime = false;

			for(int i = 0 ; i < 3 ; ++i)
			{
				UpdateLaplaciansRotation(&mesh, A, vertexDescriptors);
				SolveLinearSystem(vertexDescriptors);
			}
		}
	}
	if(g_iterateManyTimes)
	{
		for(int i = 0 ; i < 40 ; ++i)
		{
			UpdateLaplaciansRotation(&mesh, A, vertexDescriptors);
			SolveLinearSystem(vertexDescriptors);
		}
		g_iterateManyTimes = false;
	}

	if (GLpick::g_pickActive) glLoadName((GLuint)-1);

	double alpha = 1.0;

	//glEnable (GL_BLEND);
	//glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//alpha = 0.5;

	//Mesh-Faces Rendering
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL /*GL_LINE GL_FILL GL_POINT*/);
	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1., 1.);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable( GL_COLOR_MATERIAL );
	if (GLpick::g_pickActive) glLoadName((GLuint)10);

	glColor4d( 1, 1, 1, alpha );
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, &vertexDescriptors.deformedPosition[0]);
	glNormalPointer(GL_DOUBLE, 0, &vertexDescriptors.normal[0]);
	// draw the model
	glDrawElements(GL_TRIANGLES, triangles.get_size(), GL_UNSIGNED_INT, &triangles.faces[0]);
	// deactivate vertex arrays after drawing
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	if(g_showWires)
	{
		if (!GLpick::g_pickActive)
		{
			//Mesh-Edges Rendering (superimposed to faces)
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE /*GL_LINE GL_FILL GL_POINT*/);
			glColor4d( .5, .5, .5, alpha );
			glDisable( GL_LIGHTING );
				glEnableClientState(GL_VERTEX_ARRAY);
				glVertexPointer(3, GL_DOUBLE, 0, &vertexDescriptors.deformedPosition[0]);
				// draw the model
				glDrawElements(GL_TRIANGLES, triangles.get_size(), GL_UNSIGNED_INT, &triangles.faces[0]);
				// deactivate vertex arrays after drawing
				glDisableClientState(GL_VERTEX_ARRAY);
			glEnable( GL_LIGHTING );
		}
	}
	glDisable( GL_COLOR_MATERIAL );
	glDisable(GL_POLYGON_OFFSET_FILL);

	//glDisable (GL_BLEND);

	if(g_showSpheres)
	{
		//Handles rendering
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL /*GL_LINE GL_FILL GL_POINT*/);

		float	turcoise[] = { .0f, .5f, .5f, 0.3f };
		float	red[] = { .5f, .0f, .0f, 0.3f };

		for( int k = 0 ; k < handles.size() ; ++k)
		{
			if (GLpick::g_pickActive) glLoadName((GLuint)k);
			TRversor R = handles[k]->GetTRVersor();
			DrawTransparentDualSphere( _dualSphere( R * handles[k]->dS * inverse(R) ), turcoise );
		}		
	}

	glPopMatrix();

	glutSwapBuffers();
}

void reshape(GLint width, GLint height)
{
	g_viewportWidth = width;
	g_viewportHeight = height;

	// redraw viewport
	glutPostRedisplay();
}

vectorE3GA mousePosToVector(int x, int y) {
	x -= g_viewportWidth / 2;
	y -= g_viewportHeight / 2;
	return _vectorE3GA((float)-x * e1 - (float)y * e2);
}

void MouseButton(int button, int state, int x, int y)
{
	g_rotateModel = false;
	g_rotateKeyRotors = false;
	g_translateKeyRotors = false;

	if (button == GLUT_LEFT_BUTTON)
	{
		g_prevMousePos = mousePosToVector(x, y);

		GLpick::g_pickWinSize = 1;
		g_dragObject = pick(x, g_viewportHeight - y, display, &g_dragDistance);

		if(g_dragObject == -1 || g_dragObject == 10 )
		{
			vectorE3GA mousePos = mousePosToVector(x, y);
			g_rotateModel = true;

			if ((_Float(norm_e(mousePos)) / _Float(norm_e(g_viewportWidth * e1 + g_viewportHeight * e2))) < 0.2)
				g_rotateModelOutOfPlane = true;
			else g_rotateModelOutOfPlane = false;
		}
		else if(g_dragObject >= 0 && g_dragObject < handles.size())
		{
			g_rotateKeyRotors = true;
		}
	}

	if (button == GLUT_RIGHT_BUTTON)
	{
		g_prevMousePos = mousePosToVector(x, y);

		GLpick::g_pickWinSize = 1;
		g_dragObject = pick(x, g_viewportHeight - y, display, &g_dragDistance);

		if(g_dragObject >= 0 && g_dragObject < handles.size())
			g_translateKeyRotors = true;
	}
}

void MouseMotion(int x, int y)
{
	if (g_rotateModel || g_rotateKeyRotors || g_translateKeyRotors )
	{
		// get mouse position, motion
		vectorE3GA mousePos = mousePosToVector(x, y);
		vectorE3GA motion = mousePos - g_prevMousePos;

		if (g_rotateModel)
		{
			// update rotor
			if (g_rotateModelOutOfPlane)
				g_modelRotor = exp(g_camera.rotateVel * (motion ^ e3) ) * g_modelRotor;
			else 
				g_modelRotor = exp(0.00001f * (motion ^ mousePos) ) * g_modelRotor;
		}
		if(g_rotateKeyRotors)
		{
			//rotor R1 =  _rotor( inverse(g_modelRotor) * exp(-g_camera.rotateVel * (motion ^ e3) ) * g_modelRotor);
			rotor R1 =  _rotor( exp(-g_camera.rotateVel * (motion ^ e3) ) );
			if(g_dragObject < handles.size())
			{
				TRversor R = handles[g_dragObject]->R;
				handles[g_dragObject]->R = normalize(_TRversor( R1 * R  ) );
			}
		}

		if(g_translateKeyRotors)
		{
			normalizedTranslator T1 = _normalizedTranslator(inverse(g_modelRotor) * exp( _freeVector(-g_camera.translateVel*motion*ni) ) * g_modelRotor);
			if(g_dragObject < handles.size())
			{
				TRversor T = handles[g_dragObject]->T;
				handles[g_dragObject]->T = normalize(_TRversor( T1 * T ));
			}
		}

		// remember mouse pos for next motion:
		g_prevMousePos = mousePos;

		// redraw viewport
		glutPostRedisplay();
	}
}

void SpecialFunc(int key, int x, int y)
{
	switch(key) {
		case GLUT_KEY_F1 :
			{
				int mod = glutGetModifiers();
				if(mod == GLUT_ACTIVE_CTRL || mod == GLUT_ACTIVE_SHIFT )
				{
				}
			}
			break;
		case GLUT_KEY_UP:
			{
				if(g_rotateKeyRotors)
				{
					handles[g_dragObject]->dS = ChangeDualSphereRadiusSize(handles[g_dragObject]->dS, 0.025);

					// redraw viewport
					glutPostRedisplay();
				}
			}
			break;
		case GLUT_KEY_DOWN:
			{
				if(g_rotateKeyRotors)
				{
					handles[g_dragObject]->dS = ChangeDualSphereRadiusSize(handles[g_dragObject]->dS, -0.025);

					// redraw viewport
					glutPostRedisplay();
				}
			}
			break;
	}
}


void SpecialUpFunc(int key, int x, int y)
{
}

void KeyboardUpFunc(unsigned char key, int x, int y)
{
	if(key == 'w' || key == 'W')
	{
		g_showWires = !g_showWires;
		glutPostRedisplay();
	}
	
	if( key == 'h' || key == 'H' )
	{
		g_showSpheres = !g_showSpheres;
		glutPostRedisplay();
	}

	if( key == 'x' || key == 'X' )
	{
		g_iterateManyTimes = true;
		glutPostRedisplay();
	}
}

void Idle()
{
	// redraw viewport
	//glutPostRedisplay();
}

void DestroyWindow()
{
	ReleaseDrawing();
}

