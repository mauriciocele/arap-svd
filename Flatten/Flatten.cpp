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
	std::vector<Eigen::Vector3d> deformedPositions; //deformed mesh positions
	std::map<int,Eigen::Vector3d> constrainedPositions; //positional constraints
	std::vector<Eigen::Vector3d> laplacianCoordinates; //laplacian Coordinates
	std::vector<Eigen::Vector3d> normals; //for rendering (lighting)
	std::vector<Eigen::Vector3d> normalsOrig; //original normals
	std::vector<Eigen::Quaterniond> rotors;
	int size;

	VertexBuffer() : size(0)
	{
	}

	void resize(int size)
	{
		this->size = size;
		deformedPositions.resize(size);
		laplacianCoordinates.resize(size);
		normals.resize(size);
		normalsOrig.resize(size);
		rotors.resize(size);
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
	rotor R;
	translator T;
	translator Tcenter;
	dualSphere dS;
	std::set<int> constraints;

	Handle() {
		T = _rotor(1.0);
		R = _translator(1.0);
		Tcenter = _translator(1.0);
	}
	Handle(dualSphere dS)
	{
		normalizedPoint x = DualSphereCenter(dS);
		Tcenter = exp( -0.5*_vectorE3GA(x)*ni );
		T = _translator(1.0);
		R = _rotor(1.0);
		this->dS = dS;
	}

	TRversor GetTRVersor()
	{
		return _TRversor(T * _TRversor(Tcenter * R * inverse(Tcenter)));
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
std::vector<Handle> handles;


Eigen::Affine3d MotorToMatrix(const TRversor &R) {
	TRversor Ri = inverse(R);

	// compute images of basis vectors:
	c3ga::flatPoint imageOfE1NI = _flatPoint(R * c3ga::e1ni * Ri);
	c3ga::flatPoint imageOfE2NI = _flatPoint(R * c3ga::e2ni * Ri);
	c3ga::flatPoint imageOfE3NI = _flatPoint(R * c3ga::e3ni * Ri);
	c3ga::flatPoint imageOfNONI = _flatPoint(R * c3ga::noni * Ri);

	// create matrix representation:
	Eigen::Affine3d M;
	M(0, 0) = imageOfE1NI.m_c[0];
	M(1, 0) = imageOfE1NI.m_c[1];
	M(2, 0) = imageOfE1NI.m_c[2];
	M(3, 0) = imageOfE1NI.m_c[3];
	M(0, 1) = imageOfE2NI.m_c[0];
	M(1, 1) = imageOfE2NI.m_c[1];
	M(2, 1) = imageOfE2NI.m_c[2];
	M(3, 1) = imageOfE2NI.m_c[3];
	M(0, 2) = imageOfE3NI.m_c[0];
	M(1, 2) = imageOfE3NI.m_c[1];
	M(2, 2) = imageOfE3NI.m_c[2];
	M(3, 2) = imageOfE3NI.m_c[3];
	M(0, 3) = imageOfNONI.m_c[0];
	M(1, 3) = imageOfNONI.m_c[1];
	M(2, 3) = imageOfNONI.m_c[2];
	M(3, 3) = imageOfNONI.m_c[3];
	return M;
}


void ComputeLaplacianCoordinates(std::shared_ptr<SparseMatrix> A, Mesh* mesh, std::vector<Eigen::Vector3d>& laplacianCoordinates)
{
	std::fill(laplacianCoordinates.begin(), laplacianCoordinates.end(), Eigen::Vector3d(0,0,0));

	auto numRows = A->numRows();

	for( int i = 0; i < numRows ; ++i)
	{
		SparseMatrix::RowIterator aIter = A->iterator(i);
		for( ; !aIter.end() ; ++aIter )
		{
			auto j = aIter.columnIndex();
			laplacianCoordinates[i] += mesh->vertexAt(j).p * aIter.value();
		}
	}
}

bool is_constrained(std::vector<Handle>& handles, int vertex)
{
	for(Handle& handle : handles) {
		if(handle.constraints.find(vertex) != handle.constraints.end())
			return true;
	}
	return false;
}

void PreFactor(std::shared_ptr<SparseMatrix> A, std::vector<Handle>& handles)
{
	Lc = Eigen::SparseMatrix<double>(A->numRows(), A->numColumns());

	auto numRows = A->numRows();
	for( int i = 0; i < numRows ; ++i)
	{
		if(!is_constrained(handles, i))
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

	//cactus.obj
	handles.push_back(Handle(_dualSphere(c3gaPoint(.0, .04, 0.7) - 0.5*SQR(0.15)*ni)));
	handles.push_back(Handle(_dualSphere(c3gaPoint(.0, .04, -0.8) - 0.5*SQR(0.15)*ni)));

	////cylinder.obj
	//handles.push_back( boost::shared_ptr<Handle>( new Handle( _dualSphere(c3gaPoint(.0, 0.9,.0) - 0.5*SQR(0.25)*ni), false, P1 ) ) );
	//handles.push_back( boost::shared_ptr<Handle>( new Handle( _dualSphere(c3gaPoint(.0,-0.9,.0) - 0.5*SQR(0.25)*ni), true, P2 ) ) );

	//Armadillo Pie y Mano
	//handles.push_back( boost::shared_ptr<Handle>( new Handle( _dualSphere(c3gaPoint(-.5, 0.45,-.3) - 0.5*SQR(0.15)*ni), false, P1 ) ) );
	//handles.push_back( boost::shared_ptr<Handle>( new Handle( _dualSphere(c3gaPoint(.2,-0.6,.1) - 0.5*SQR(0.15)*ni), true, P2 ) ) );
	//Armadillo Pubis y Cabeza
	//handles.push_back( boost::shared_ptr<Handle>( new Handle( _dualSphere(c3gaPoint(.0, 0.4,-.2) - 0.5*SQR(0.15)*ni), false, P1 ) ) );
	//handles.push_back( boost::shared_ptr<Handle>( new Handle( _dualSphere(c3gaPoint(.0,-0.05,.1) - 0.5*SQR(0.15)*ni), true, P2 ) ) );

	vertexDescriptors.resize(mesh.numVertices());
	triangles.resize(mesh.numFaces()*3);
	int n = vertexDescriptors.get_size();

	for(Vertex& vertex : mesh.getVertices()) {
		vertexDescriptors.normalsOrig[vertex.ID] = vertex.n;

		normalizedPoint position = c3gaPoint( vertex.p.x(), vertex.p.y(), vertex.p.z() );

		for( Handle& handle : handles)
		{
			TRversor TR = handle.GetTRVersor();

			if( _double(position << (TR * handle.dS * inverse(TR))) > 0 ) //inside the sphere
			{
				handle.constraints.insert(vertex.ID);
				vertexDescriptors.constrainedPositions[vertex.ID] = vertex.p;
			}
		}
	}

	for(Face& face : mesh.getFaces()) {
		int i = face.ID;
		int	v1 = face.edge->vertex->ID;
		int	v2 = face.edge->next->vertex->ID;
		int	v3 = face.edge->next->next->vertex->ID;
		triangles.faces[i*3 + 0] = v1;
		triangles.faces[i*3 + 1] = v2;
		triangles.faces[i*3 + 2] = v3;
	}

	A = CreateLaplacianMatrix( &mesh, systemType );
	
	ComputeLaplacianCoordinates(A, &mesh, vertexDescriptors.laplacianCoordinates);

	b3 = Eigen::MatrixXd(A->numRows(), 3);

	PreFactor(A, handles);

	glutMainLoop();

	return 0;
}

void SolveLinearSystem(VertexBuffer& vertexDescriptors)
{
	int n = vertexDescriptors.get_size();
	for( int i = 0 ; i < n ; ++i ) {
		b3.row(i) = vertexDescriptors.laplacianCoordinates[i];
	}

	for( Handle& handle : handles) {
		Eigen::Affine3d M = MotorToMatrix(handle.GetTRVersor());
		for(int i : handle.constraints) {
			b3.row(i) = M * vertexDescriptors.constrainedPositions[i];
		}
	}

	xyz = solver.solve(b3);

	for( int i = 0 ; i < n ; ++i )
	{
		vertexDescriptors.deformedPositions[i] = xyz.row(i);
		vertexDescriptors.normals[i] = vertexDescriptors.rotors[i]._transformVector(vertexDescriptors.normalsOrig[i]);
	}
}

void UpdateLaplaciansRotation(Mesh *mesh, std::shared_ptr<SparseMatrix> A, VertexBuffer& vertexDescriptors)
{
	Eigen::Matrix3d m;
	for( Vertex& vertex : mesh->getVertices() )
	{
		m.setZero();
		int i = vertex.ID;
		const Eigen::Vector3d &pi = mesh->vertexAt(i).p;
		const Eigen::Vector3d &tpi = vertexDescriptors.deformedPositions[i];
		double S = 0;
		for(Vertex::EdgeAroundIterator edgeAroundIter = vertex.iterator() ; !edgeAroundIter.end() ; edgeAroundIter++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;
			const Eigen::Vector3d &pj = mesh->vertexAt(j).p;
			const Eigen::Vector3d &tpj = vertexDescriptors.deformedPositions[j];
			Eigen::Vector3d eij = (pj - pi) * (*A)(i, j);
			Eigen::Vector3d teij = tpj - tpi;
			m += eij * teij.transpose();
			S += eij.dot(eij);
			S += teij.dot(teij) * (*A)(i, j);
		}
		vertexDescriptors.rotors[i] = GARotorEstimator(m, S);
	}

	std::fill(vertexDescriptors.laplacianCoordinates.begin(), vertexDescriptors.laplacianCoordinates.end(), Eigen::Vector3d(0,0,0));
	for(Vertex& vertex : mesh->getVertices())
	{
		int i = vertex.ID;
		for(Vertex::EdgeAroundIterator edgeAroundIter = vertex.iterator() ; !edgeAroundIter.end() ; edgeAroundIter++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;
			double wij = (*A)(i, j);
			Eigen::Quaterniond &Ri = vertexDescriptors.rotors[i];
			Eigen::Quaterniond &Rj = vertexDescriptors.rotors[j];
			Eigen::Vector3d V = mesh->vertexAt(j).p - mesh->vertexAt(i).p;
			vertexDescriptors.laplacianCoordinates[i] += 0.5 * wij * (Ri._transformVector(V) + Rj._transformVector(V));
		}
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

			for(int i = 0 ; i < 10 ; ++i)
			{
				UpdateLaplaciansRotation(&mesh, A, vertexDescriptors);
				SolveLinearSystem(vertexDescriptors);
			}
		}
	}
	if(g_iterateManyTimes)
	{
		for(int i = 0 ; i < 400 ; ++i)
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
	glVertexPointer(3, GL_DOUBLE, 0, &vertexDescriptors.deformedPositions[0]);
	glNormalPointer(GL_DOUBLE, 0, &vertexDescriptors.normals[0]);
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
				glVertexPointer(3, GL_DOUBLE, 0, &vertexDescriptors.deformedPositions[0]);
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
			TRversor R = handles[k].GetTRVersor();
			DrawTransparentDualSphere( _dualSphere( R * handles[k].dS * inverse(R) ), turcoise );
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
				rotor R = handles[g_dragObject].R;
				handles[g_dragObject].R = normalize(_TRversor( R1 * R  ) );
			}
		}

		if(g_translateKeyRotors)
		{
			normalizedTranslator T1 = _normalizedTranslator(inverse(g_modelRotor) * exp( _freeVector(-g_camera.translateVel*motion*ni) ) * g_modelRotor);
			if(g_dragObject < handles.size())
			{
				translator T = handles[g_dragObject].T;
				handles[g_dragObject].T = normalize(_TRversor( T1 * T ));
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
					handles[g_dragObject].dS = ChangeDualSphereRadiusSize(handles[g_dragObject].dS, 0.025);

					// redraw viewport
					glutPostRedisplay();
				}
			}
			break;
		case GLUT_KEY_DOWN:
			{
				if(g_rotateKeyRotors)
				{
					handles[g_dragObject].dS = ChangeDualSphereRadiusSize(handles[g_dragObject].dS, -0.025);

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

