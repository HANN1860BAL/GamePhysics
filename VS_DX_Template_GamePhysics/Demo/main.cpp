//--------------------------------------------------------------------------------------
// File: main.cpp
//
// The main file containing the entry point main().
//--------------------------------------------------------------------------------------

#include <sstream>
#include <iomanip>
#include <random>
#include <iostream>

//DirectX includes
#include <DirectXMath.h>
using namespace DirectX;

// Effect framework includes
#include <d3dx11effect.h>

// DXUT includes
#include <DXUT.h>
#include <DXUTcamera.h>

// DirectXTK includes
#include "Effects.h"
#include "VertexTypes.h"
#include "PrimitiveBatch.h"
#include "GeometricPrimitive.h"
#include "ScreenGrab.h"

// AntTweakBar includes
#include "AntTweakBar.h"

// Internal includes
#include "util/util.h"
#include "util/FFmpeg.h"

#include <vector>
#include "util/collisionDetect.h"
#include <string>

// DXUT camera
// NOTE: CModelViewerCamera does not only manage the standard view transformation/camera position 
//       (CModelViewerCamera::GetViewMatrix()), but also allows for model rotation
//       (CModelViewerCamera::GetWorldMatrix()). 
//       Look out for CModelViewerCamera::SetButtonMasks(...).
CModelViewerCamera g_camera;

// Effect corresponding to "effect.fx"
ID3DX11Effect* g_pEffect = nullptr;

// Main tweak bar
TwBar* g_pTweakBar;

// DirectXTK effects, input layouts and primitive batches for different vertex types
BasicEffect*                               g_pEffectPositionColor = nullptr;
ID3D11InputLayout*                         g_pInputLayoutPositionColor = nullptr;
PrimitiveBatch<VertexPositionColor>*       g_pPrimitiveBatchPositionColor = nullptr;

BasicEffect*                               g_pEffectPositionNormal = nullptr;
ID3D11InputLayout*                         g_pInputLayoutPositionNormal = nullptr;
PrimitiveBatch<VertexPositionNormal>*      g_pPrimitiveBatchPositionNormal = nullptr;

BasicEffect*                               g_pEffectPositionNormalColor = nullptr;
ID3D11InputLayout*                         g_pInputLayoutPositionNormalColor = nullptr;
PrimitiveBatch<VertexPositionNormalColor>* g_pPrimitiveBatchPositionNormalColor = nullptr;

// DirectXTK simple geometric primitives
std::unique_ptr<GeometricPrimitive> g_pSphere;
std::unique_ptr<GeometricPrimitive> g_pTeapot;
std::unique_ptr<GeometricPrimitive> g_pCube;

// Movable object management
XMINT2   g_viMouseDelta = XMINT2(0, 0);
XMFLOAT3 g_vfMovableObjectPos = XMFLOAT3(0, 0, 0);

// TweakAntBar GUI variables

//Multiple usage
int   g_iNumSpheres = 3; //number of spheres in a row
float g_fSphereSize = 0.05f;
float g_fDamping = 0.001f;
float g_fMass = 1.0f;

//MassSpringSystem
bool g_bMassSpringSystem = false;
bool g_bEuler = false;
bool g_bMidpoint = false;
bool g_bRungeKutta = false;
bool g_bFloorCollsion = false;
bool g_bSphereCollsion = false;
float g_fStiffness = 10.0f;
float g_fTimeStepSize = 0.01f;

//Rigidbody simulation
bool g_bRigidbody = false;
bool g_bInteractionLeft = false;
bool g_bInteractionRight = false;
float g_fImpulseConstant = 0.5f;

//Balls in a box
bool g_bBiab = false;
bool g_bBiabNaive = false;
bool g_bBiabKDTree = false;
bool g_bBiabUniformGrid = false;
float g_fcollisionScalar = 25.0f;

//Global
bool b_start = true;
float f_gravity = -9.81f;
float f_oldSize = 0;
static double f_timeAcc;
int i_oldNum = 0;
ID3D11Device* pd3dDeviceTweakbar;

//usefull for initialisation
XMVECTOR XMV_zero = XMVectorSet(0.0f, 0.0f, 0.0f, 0.0f);
XMMATRIX XMM_zero = XMMatrixSet(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
XMMATRIX XMM_identity = XMMatrixSet(1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f);

// Video recorder
FFmpeg* g_pFFmpegVideoRecorder = nullptr;

//Points for MassSpringSystem
struct Point
{
	XMVECTOR XMV_position;
	XMVECTOR XMV_velocity;
	XMVECTOR XMV_force;
	float f_mass;
	bool b_Static;
	bool b_Dummy;
	int i_id;
};

//Springs fpr MassSpringSystem
struct Spring
{
	int i_point1;
	int i_point2;;
	float f_stiffness;
	float f_initLength;
	float f_currentLength;
};

//corners of Rigidbodys
struct Corner
{
	XMVECTOR XMV_position;
	XMVECTOR XMV_velocity;
	XMVECTOR XMV_force;
};

//rigidbody boxes
struct Box
{
	XMVECTOR XMV_position;
	XMMATRIX XMM_inertiaTensor;
	XMVECTOR XMV_orientation;
	XMVECTOR XMV_linearVelocity;
	XMVECTOR XMV_angularVelocity;
	XMVECTOR XMV_angularMomentum;
	XMVECTOR XMV_forceAccumulator;
	XMVECTOR XMV_torqueAccumulator;
	XMMATRIX XMM_transform;
	std::vector<Corner> v_corner;
	float f_mass;
	float f_lengthX;
	float f_lengthY;
	float f_lengthZ;
	Box(XMVECTOR XMV_position, float f_mass, float f_lengthX, float f_lengthY, float f_lengthZ)
		: XMV_position(XMV_position), f_mass(f_mass), f_lengthX(f_lengthX), f_lengthY(f_lengthY), f_lengthZ(f_lengthZ){}
};

//Balls for MassCollisionDetection
struct Ball
{
	XMVECTOR XMV_position;
	XMVECTOR XMV_velocity;
	int id;
};

//Helpstruct for MassCollisiondetection
struct Pairs
{
	Ball* ba_ballA;
	Ball* ba_ballB;
};

//storage for MassSpringSystem
std::vector<Point> v_point;
std::vector<Spring> v_spring;

//storage for Rigidbodys
std::vector<Box> v_box;

//storage for MassCollisionDetection
std::vector<Ball> v_ball;
std::vector<std::vector<Ball>> vv_kdtree;
std::vector<Pairs> v_pairs;


//Utility
Point GetPointOf(int id)
{
	for (int i = 0; i < v_point.size(); i++)
	{
		if (id == v_point[i].i_id)
			return v_point[i];
	}
}

//shows XMMATRIX on the console
void PrintMatrix(XMMATRIX XMM_matrix, std::string s_name)
{
	XMFLOAT4X4 XMF4X4_tmp;
	XMStoreFloat4x4(&XMF4X4_tmp, XMM_matrix);
	std::cout << s_name << "\n";
	std::cout << std::setw(15) << XMF4X4_tmp._11 << " " << std::setw(15) << XMF4X4_tmp._12 << " " << std::setw(15) << XMF4X4_tmp._13 << " " << std::setw(15) << XMF4X4_tmp._14 << "\n";
	std::cout << std::setw(15) << XMF4X4_tmp._21 << " " << std::setw(15) << XMF4X4_tmp._22 << " " << std::setw(15) << XMF4X4_tmp._23 << " " << std::setw(15) << XMF4X4_tmp._24 << "\n";
	std::cout << std::setw(15) << XMF4X4_tmp._31 << " " << std::setw(15) << XMF4X4_tmp._32 << " " << std::setw(15) << XMF4X4_tmp._33 << " " << std::setw(15) << XMF4X4_tmp._34 << "\n";
	std::cout << std::setw(15) << XMF4X4_tmp._41 << " " << std::setw(15) << XMF4X4_tmp._42 << " " << std::setw(15) << XMF4X4_tmp._43 << " " << std::setw(15) << XMF4X4_tmp._44 << "\n";
}

//shows XMVECTOR on the console
void PrintVector(XMVECTOR XMV_vector, std::string s_name)
{
	XMFLOAT4 XMF4_tmp;
	XMStoreFloat4(&XMF4_tmp, XMV_vector);
	std::cout << s_name << "\n";
	std::cout << "w: " << std::setw(15) << XMF4_tmp.w << " x: " << std::setw(15) << XMF4_tmp.x << " y: " << std::setw(15) << XMF4_tmp.y << " z: " << std::setw(15) << XMF4_tmp.z << "\n";
}

//transforms a XMVECTOR(w,x,y,z) to a quaternion(x,y,z,w)
XMVECTOR VectorToQuaternion(XMVECTOR XMV_vector)
{
	XMFLOAT4 XMF4_tmp;
	XMStoreFloat4(&XMF4_tmp, XMV_vector);
	XMV_vector = XMVectorSet(XMF4_tmp.x, XMF4_tmp.y, XMF4_tmp.z, XMF4_tmp.w);
	return XMV_vector;
}

//transposes a XMVECTOR and multiplys with itself
XMMATRIX VectorMulTranspose(XMVECTOR XMV_vector)
{
	std::vector<float> v_Ovector;
	std::vector<float> v_Tvector;
	std::vector<float> v_Mmatrix;

	v_Ovector.push_back(XMVectorGetX(XMV_vector));
	v_Ovector.push_back(XMVectorGetY(XMV_vector));
	v_Ovector.push_back(XMVectorGetZ(XMV_vector));
	v_Ovector.push_back(XMVectorGetW(XMV_vector));

	v_Tvector.push_back(XMVectorGetX(XMV_vector));
	v_Tvector.push_back(XMVectorGetY(XMV_vector));
	v_Tvector.push_back(XMVectorGetZ(XMV_vector));
	v_Tvector.push_back(XMVectorGetW(XMV_vector));

	for (int i = 0; i < v_Ovector.size(); i++)
	{
		for (int k = 0; k < v_Tvector.size(); k++)
		{
			v_Mmatrix.push_back(v_Ovector[i] * v_Tvector[k]);
		}
	}

	return XMMatrixSet(v_Mmatrix[0], v_Mmatrix[1], v_Mmatrix[2], v_Mmatrix[3],
		v_Mmatrix[4], v_Mmatrix[5], v_Mmatrix[6], v_Mmatrix[7],
		v_Mmatrix[8], v_Mmatrix[9], v_Mmatrix[10], v_Mmatrix[11],
		v_Mmatrix[12], v_Mmatrix[13], v_Mmatrix[14], v_Mmatrix[15]);
}

//Resets the simulation
void Reset()
{
	b_start = true;
	g_bEuler = false;
	g_bMidpoint = false;
	g_bRungeKutta = false;
	std::cout << "-----------------------------------------\n";
}

//Setter and functions for the twearkbaroptions
void SetEuler()
{
	g_bEuler = true;
	g_bMidpoint = false;
	g_bRungeKutta = false;
}

void SetMidpoint()
{
	g_bMidpoint = true;
	g_bEuler = false;
	g_bRungeKutta = false;
}

void SetRungeKutta()
{
	g_bRungeKutta = true;
	g_bMidpoint = false;
	g_bEuler = false;
}

void SetMassSpringSystem()
{
	b_start = true;

	g_bMassSpringSystem = true;
	g_bRigidbody = false;
	g_bBiab = false;
	v_box.clear();
	v_ball.clear();
}

void SetRigidBody()
{
	b_start = true;

	g_bRigidbody = true;
	g_bMassSpringSystem = false;
	g_bBiab = false;
	v_point.clear();
	v_spring.clear();
	v_ball.clear();
}

void SetBiab()
{
	b_start = true;

	g_bBiab = true;
	g_bMassSpringSystem = false;
	g_bRigidbody = false;
	v_box.clear();
	v_point.clear();
	v_spring.clear();	
}

void SetBiabNaive()
{
	g_bBiabNaive = true;
	g_bBiabKDTree = false;
	g_bBiabUniformGrid = false;
}

void SetBiabKDTree()
{
	g_bBiabNaive = false;
	g_bBiabKDTree = true;
	g_bBiabUniformGrid = false;
}

void SetBiabUniformGrid()
{
	g_bBiabNaive = false;
	g_bBiabKDTree = false;
	g_bBiabUniformGrid = true;
}

void SetStiffness(float f_stiffness)
{
	if (v_spring[0].f_stiffness != f_stiffness){
		for (int i = 0; i < v_spring.size(); i++)
		{
			v_spring[i].f_stiffness = f_stiffness;
		}
	}
}

void SetMass(float f_mass)
{
	if (v_point[0].f_mass != f_mass)
	{
		for (int i = 0; i < v_point.size(); i++)
		{
			v_point[i].f_mass = f_mass;
		}
	}
}

float RandomBetween(float smallNumber, float bigNumber)
{
	float diff = bigNumber - smallNumber;
	return (((float)rand() / RAND_MAX) * diff) + smallNumber;
}

void randomPosition()
{
	int i = (int)RandomBetween(0, (int)v_point.size() - 1);
	if (!v_point[i].b_Static)
	{
		XMFLOAT3 XMF_tmp;
		XMStoreFloat3(&XMF_tmp, v_point[i].XMV_position);
		XMF_tmp.y += RandomBetween(-0.5f, 0.5f);
		v_point[i].XMV_position = XMLoadFloat3(&XMF_tmp);
	}
	else randomPosition();
}

void InteractionLeft()
{
	XMFLOAT3 XMF3_tmp = XMFLOAT3(-0.1f, 0.0f, 0.0f);
	v_box[1].XMV_position += XMLoadFloat3(&XMF3_tmp);
}

void InteractionRight()
{
	XMFLOAT3 XMF3_tmp = XMFLOAT3(0.1f, 0.0f, 0.0f);
	v_box[1].XMV_position += XMLoadFloat3(&XMF3_tmp);
}


// Create TweakBar and add required buttons and variables
void InitTweakBar(ID3D11Device* pd3dDevice)
{
	TwInit(TW_DIRECT3D11, pd3dDevice);

	g_pTweakBar = TwNewBar("TweakBar");

	// HINT: For buttons you can directly pass the callback function as a lambda expression.
	TwAddButton(g_pTweakBar, "Reset Camera", [](void *){g_camera.Reset(); }, nullptr, "");
	TwAddButton(g_pTweakBar, "Mass Spring System", [](void *){SetMassSpringSystem(); }, nullptr, "");
	TwAddButton(g_pTweakBar, "Rigidbody Simulation", [](void *){SetRigidBody(); }, nullptr, "");
	TwAddButton(g_pTweakBar, "Balls in a box", [](void *){SetBiab(); }, nullptr, "");
	TwAddButton(g_pTweakBar, "Reset", [](void *){Reset(); }, nullptr, "");
	TwAddVarRW(g_pTweakBar, "Time Step Size", TW_TYPE_FLOAT, &g_fTimeStepSize, "min=0.001 step=0.001");
	TwAddSeparator(g_pTweakBar, "------------------------", "");

	if (g_bMassSpringSystem)
	{
		TwAddButton(g_pTweakBar, "Euler", [](void *){SetEuler(); }, nullptr, "");
		TwAddButton(g_pTweakBar, "Midpoint", [](void *){SetMidpoint(); }, nullptr, "");
		TwAddButton(g_pTweakBar, "RungeKutta", [](void *){SetRungeKutta(); }, nullptr, "");
		TwAddButton(g_pTweakBar, "Random position", [](void *){randomPosition(); }, nullptr, "");
		TwAddVarRW(g_pTweakBar, "Floor Collsions", TW_TYPE_BOOLCPP, &g_bFloorCollsion, "");
		TwAddVarRW(g_pTweakBar, "Sphere Collsions", TW_TYPE_BOOLCPP, &g_bSphereCollsion, "");
		TwAddVarRW(g_pTweakBar, "Num Spheres", TW_TYPE_INT32, &g_iNumSpheres, "min=1");
		TwAddVarRW(g_pTweakBar, "Sphere Size", TW_TYPE_FLOAT, &g_fSphereSize, "min=0.01 step=0.01");
		TwAddVarRW(g_pTweakBar, "Sphere Mass", TW_TYPE_FLOAT, &g_fMass, "min=0.01 step=0.01");
		TwAddVarRW(g_pTweakBar, "Spring Stiffness", TW_TYPE_FLOAT, &g_fStiffness, " ");
		TwAddVarRW(g_pTweakBar, "Damping", TW_TYPE_FLOAT, &g_fDamping, "max=0 step=0.01");
	}

	if (g_bRigidbody)
	{
		TwAddButton(g_pTweakBar, "Interaction left", [](void*){InteractionLeft(); }, nullptr, "");
		TwAddButton(g_pTweakBar, "Interaction right", [](void*){InteractionRight(); }, nullptr, "");
	}

	if (g_bBiab)
	{
		TwAddButton(g_pTweakBar, "Biab: Naive", [](void *){SetBiabNaive(); }, nullptr, "");
		TwAddButton(g_pTweakBar, "Biab: KD Tree", [](void *){SetBiabKDTree(); }, nullptr, "");
		TwAddButton(g_pTweakBar, "Biab: Uniform grid", [](void *){SetBiabUniformGrid(); }, nullptr, "");
		TwAddVarRW(g_pTweakBar, "Num Spheres", TW_TYPE_INT32, &g_iNumSpheres, "min=1");
		TwAddVarRW(g_pTweakBar, "Sphere Size", TW_TYPE_FLOAT, &g_fSphereSize, "min=0.01 step=0.01");
		TwAddVarRW(g_pTweakBar, "Sphere Mass", TW_TYPE_FLOAT, &g_fMass, "min=0.001 step=0.01");
		TwAddVarRW(g_pTweakBar, "Collision Scalar", TW_TYPE_FLOAT, &g_fcollisionScalar, "min=0 step=1");
		TwAddVarRW(g_pTweakBar, "Damping", TW_TYPE_FLOAT, &g_fDamping, "min=0 max=1");
	}
}



//////////////MASSSPRINGSYSTEM/////////////////////
void InitPoints()
{
	v_point.clear();

	float f_distance = 1.0f / (g_iNumSpheres - 1);
	float f_offSet = -0.5f;
	if (f_distance < g_fSphereSize * 2)
	{
		f_distance = g_fSphereSize * 2;
		f_offSet -= f_distance;
	}
	int i_count = 0;
	for (int i = 0; i < g_iNumSpheres + 2; i++){
		for (int k = 0; k < g_iNumSpheres + 2; k++){

			Point p_point = Point();
			p_point.b_Dummy = false;
			p_point.b_Static = false;
			XMFLOAT3 XMF3_position = XMFLOAT3((f_offSet - f_distance) + k*f_distance, 0.5f, (f_offSet - f_distance) + i*f_distance);

			p_point.XMV_position = XMLoadFloat3(&XMF3_position);
			p_point.XMV_force = XMLoadFloat3(&XMFLOAT3(0.0f, 0.0f, 0.0f));
			p_point.XMV_velocity = p_point.XMV_force;
			p_point.f_mass = g_fMass;

			//set staticflag
			if (k == 1 || i == 1 || k == g_iNumSpheres || i == g_iNumSpheres)
				p_point.b_Static = true;
			//set dummyflag
			if (k == 0 || i == 0 || k == g_iNumSpheres + 1 || i == g_iNumSpheres + 1)
				p_point.b_Dummy = true;

			p_point.i_id = i_count;
			i_count++;
			v_point.push_back(p_point);
		}
	}
}

void InitSprings()
{
	v_spring.clear();

	Spring sp_spring;
	sp_spring.f_stiffness = g_fStiffness;

	//set springs to all points except, but not from dummys
	for (int i = 0; i < v_point.size() - 1; i++){
		if (!v_point[i].b_Dummy){
			sp_spring.i_point1 = v_point[i].i_id;
			sp_spring.i_point2 = v_point[i + 1].i_id;
			v_spring.push_back(sp_spring);
			sp_spring.i_point2 = v_point[i - 1].i_id;
			v_spring.push_back(sp_spring);
			sp_spring.i_point2 = v_point[i + g_iNumSpheres + 1].i_id;
			v_spring.push_back(sp_spring);
			sp_spring.i_point2 = v_point[i + g_iNumSpheres + 2].i_id;
			v_spring.push_back(sp_spring);
			sp_spring.i_point2 = v_point[i + g_iNumSpheres + 3].i_id;
			v_spring.push_back(sp_spring);
			sp_spring.i_point2 = v_point[i - g_iNumSpheres - 1].i_id;
			v_spring.push_back(sp_spring);
			sp_spring.i_point2 = v_point[i - g_iNumSpheres - 2].i_id;
			v_spring.push_back(sp_spring);
			sp_spring.i_point2 = v_point[i - g_iNumSpheres - 3].i_id;
			v_spring.push_back(sp_spring);
		}
	}



	//claculates length of spring
	for (int i = 0; i < v_spring.size(); i++)
	{
		//brakets from in -> out: subtract b_start- and endposition of spring, calculate length of solution of subtraction, store solution as float in f_initLength
		XMStoreFloat(&v_spring[i].f_initLength, XMVector3Length(XMVectorSubtract(GetPointOf(v_spring[i].i_point1).XMV_position, GetPointOf(v_spring[i].i_point2).XMV_position)));
		v_spring[i].f_currentLength = v_spring[i].f_initLength;
	}

	//delete dummys and their springs
	for (int i = 0; i < v_point.size(); i++)
	{
		if (v_point[i].b_Dummy)
		{
			for (int k = 0; k < v_spring.size(); k++)
			{
				if (v_point[i].i_id == v_spring[k].i_point2)
				{
					v_spring.erase(v_spring.begin() + k);
				}
			}
			v_point.erase(v_point.begin() + i);
			i--;
		}
	}

	//delete double springs
	for (int m = 0; m < v_spring.size(); m++){
		for (int n = 0; n < v_spring.size(); n++){
			if (v_spring[m].i_point1 == v_spring[n].i_point2 && v_spring[m].i_point2 == v_spring[n].i_point1){
				v_spring.erase(v_spring.begin() + n);
			}
		}
	}
}

XMVECTOR UseEulerIntegration(Point* p_point, float f_timeStep)
{
	//apply gravity
	p_point->XMV_force = XMVectorSet(0.0f, f_gravity, 0.0f, 0.0f) * p_point->f_mass;

	//add internal forces
	for (int k = 0; k < v_spring.size(); k++)
	{
		if (p_point->i_id == v_spring[k].i_point1)
		{
			p_point->XMV_force += -v_spring[k].f_stiffness * (v_spring[k].f_currentLength - v_spring[k].f_initLength) * ((p_point->XMV_position - GetPointOf(v_spring[k].i_point2).XMV_position) / v_spring[k].f_currentLength) + g_fDamping * p_point->XMV_velocity;
		}

		if (p_point->i_id == v_spring[k].i_point2)
		{
			p_point->XMV_force += -v_spring[k].f_stiffness * (v_spring[k].f_currentLength - v_spring[k].f_initLength) * ((p_point->XMV_position - GetPointOf(v_spring[k].i_point1).XMV_position) / v_spring[k].f_currentLength) + g_fDamping * p_point->XMV_velocity;
		}
	}

	p_point->XMV_velocity += (p_point->XMV_force / p_point->f_mass)* f_timeStep; //pow(f_timeStep,2);

	//apply damping
	//p_point->XMV_velocity += g_fDamping * p_point->XMV_velocity;
	//
	////apply damping
	//p_point->XMV_velocity += g_fDamping * XMVectorAbs(p_point->XMV_velocity) * f_timeAcc;
	//
	//p_point->XMV_force *= pow(f_timeStep, 2);
	//
	//XMVECTOR XMV_newPosition = p_point->XMV_position + p_point->XMV_force;
	//
	////apply damping
	//p_point->XMV_velocity = p_point->XMV_force * f_timeStep;
	//XMStoreFloat3(&XMF3_damping, g_fDamping *p_point->XMV_velocity *f_timeAcc);
	//p_point->XMV_force += g_fDamping * p_point->XMV_velocity * f_timeAcc;
	//
	//XMFLOAT3 tmp;
	//
	//if (truncf(pow(g_iNumSpheres+2, 2) / 2.0f) == p_point->i_id)
	//{
	//	//std::cout << "gravity: " << XMF3_gravity.y << "\n";
	//	//std::cout << "internalforce: " << XMF3_internalForce.x << " " << XMF3_internalForce.y << " " << XMF3_internalForce.z << "\n";
	//	//std::cout << "damping: " << XMF3_damping.x << " " << XMF3_damping.y << " " << XMF3_damping.z << "\n";
	//	//XMStoreFloat3(&tmp, p_point->XMV_force);
	//	//std::cout << "Force: " << tmp.x << " " << tmp.y << " " << tmp.z << "\n";
	//	XMStoreFloat3(&tmp, p_point->XMV_velocity);
	//	std::cout << "Velocity: " << tmp.x << " " << tmp.y << " " << tmp.z << "\n";
	//}
	//
	//XMFLOAT3 tmp;
	//XMStoreFloat3(&tmp, p_point->XMV_velocity);
	//std::cout << "Velocity: " << tmp.x << " " << tmp.y << " " << tmp.z << "\n";
	//
	//XMStoreFloat3(&tmp, p_point->XMV_force);
	//std::cout <<"Force: " << tmp.x << " " << tmp.y << " " << tmp.z << "\n";

	return p_point->XMV_velocity;
}

void CalculateCurrentLength()
{
	for (int i = 0; i < v_spring.size(); i++)
	{
		//brakets from in -> out: subtract b_start- and endposition of spring, calculate length of solution of subtraction, store solution as float in f_initLength
		XMStoreFloat(&v_spring[i].f_currentLength, XMVector3Length(XMVectorSubtract(GetPointOf(v_spring[i].i_point1).XMV_position, GetPointOf(v_spring[i].i_point2).XMV_position)));
	}
}

//very simple Collisiondetection for MassSpringSystem
void CollisionDetectionSpheres()
{
	for (int i = 0; i < v_point.size(); i++)
	{
		for (int k = 0; k < v_point.size(); k++)
		{
			float f_tmp;
			XMStoreFloat(&f_tmp, XMVector3Length(v_point[i].XMV_position - v_point[k].XMV_position));
			if (f_tmp < 2 * g_fSphereSize)
			{
				f_tmp = 2 * g_fSphereSize - f_tmp;
				v_point[i].XMV_position += XMVector3Normalize(v_point[i].XMV_position - v_point[k].XMV_position)*(f_tmp / 2);
				v_point[k].XMV_position += XMVector3Normalize(v_point[i].XMV_position - v_point[k].XMV_position)*(-f_tmp / 2);
			}
		}
	}
}

void ApplyPhysikMSS()
{
	if (g_bEuler)
	{
		for (int i = 0; i < v_point.size(); i++)
		{
			UseEulerIntegration(&v_point[i], g_fTimeStepSize);
		}
	}
	else if (g_bMidpoint)
	{
		XMVECTOR XMV_k1;
		XMVECTOR XMV_tmp;

		for (int i = 0; i < v_point.size(); i++)
		{
			XMV_tmp = v_point[i].XMV_velocity;
			XMV_k1 = UseEulerIntegration(&v_point[i], g_fTimeStepSize / 2);
			//v_point[i].XMV_velocity = XMV_tmp + XMV_k1;
			XMV_k1 = UseEulerIntegration(&v_point[i], g_fTimeStepSize);
			v_point[i].XMV_velocity = XMV_tmp + XMV_k1;
		}
	}
	else if (g_bRungeKutta)
	{
		XMVECTOR XMV_k1;
		XMVECTOR XMV_k2;
		XMVECTOR XMV_k3;
		XMVECTOR XMV_k4;
		XMVECTOR XMV_tmp;

		for (int i = 0; i < v_point.size(); i++)
		{
			XMV_tmp = v_point[i].XMV_velocity;
			XMV_k1 = v_point[i].XMV_velocity; //UseEulerIntegration(&v_point[i], g_fTimeStepSize);
			//v_point[i].XMV_velocity = XMV_tmp + XMV_k1;
			XMV_k2 = UseEulerIntegration(&v_point[i], g_fTimeStepSize / 2);
			//v_point[i].XMV_velocity = XMV_tmp + XMV_k2;
			XMV_k3 = UseEulerIntegration(&v_point[i], g_fTimeStepSize / 2);
			//v_point[i].XMV_velocity = XMV_tmp + XMV_k3;
			XMV_k4 = UseEulerIntegration(&v_point[i], g_fTimeStepSize);
			v_point[i].XMV_velocity = XMV_tmp + (g_fTimeStepSize / 6.0f) * (XMV_k1 + 2 * XMV_k2 + 2 * XMV_k3 + XMV_k4);
		}
	}

	if (g_bEuler || g_bMidpoint || g_bRungeKutta)
	{
		if (g_bSphereCollsion)
			CollisionDetectionSpheres();
		XMFLOAT3 tmp;
		for (int i = 0; i < v_point.size(); i++)
		{
			if (!v_point[i].b_Static)
			{
				v_point[i].XMV_position += v_point[i].XMV_velocity;
				//if (truncf(pow(g_iNumSpheres + 2, 2) / 2.0f) == v_point[i].i_id)
				//{
				//	XMStoreFloat3(&tmp, v_point[i].XMV_position);
				//	std::cout << "Position y: " << tmp.y << "\n";
				//}
				if (g_bFloorCollsion)
				{
					XMStoreFloat3(&tmp, v_point[i].XMV_position);
					if (tmp.y <= -1.0f + g_fSphereSize)
						tmp.y = -1.0f + g_fSphereSize;
					v_point[i].XMV_position = XMLoadFloat3(&tmp);
				}
			}
		}
		CalculateCurrentLength();
	}
}



/////////////RIGISBODYSIMULATION///////////////////
void CalculateCorners(Box* b_box)
{
	XMFLOAT3 XMF3_tmp;
	XMStoreFloat3(&XMF3_tmp, b_box->XMV_position);
	Corner c_corner;
	c_corner.XMV_velocity = XMV_zero;
	c_corner.XMV_force = XMVectorSet(0.0f, f_gravity, 0.0f, 0.0f) * (b_box->f_mass / 8);
	//c_corner.XMV_force = XMVector4Transform(XMVectorSet(0.0f, f_gravity, 0.0f, 0.0f) * (b_box->f_mass / 8), XMMatrixInverse(NULL, XMMatrixRotationQuaternion(b_box->XMV_orientation)));
	c_corner.XMV_position = XMVectorSet(XMF3_tmp.x - b_box->f_lengthX / 2, XMF3_tmp.y - b_box->f_lengthY / 2, XMF3_tmp.z - b_box->f_lengthZ / 2, 0.0f) - b_box->XMV_position;
	b_box->v_corner.push_back(c_corner);
	c_corner.XMV_position = XMVectorSet(XMF3_tmp.x - b_box->f_lengthX / 2, XMF3_tmp.y - b_box->f_lengthY / 2, XMF3_tmp.z + b_box->f_lengthZ / 2, 0.0f) - b_box->XMV_position;
	b_box->v_corner.push_back(c_corner);
	c_corner.XMV_position = XMVectorSet(XMF3_tmp.x - b_box->f_lengthX / 2, XMF3_tmp.y + b_box->f_lengthY / 2, XMF3_tmp.z - b_box->f_lengthZ / 2, 0.0f) - b_box->XMV_position;
	b_box->v_corner.push_back(c_corner);
	c_corner.XMV_position = XMVectorSet(XMF3_tmp.x - b_box->f_lengthX / 2, XMF3_tmp.y + b_box->f_lengthY / 2, XMF3_tmp.z + b_box->f_lengthZ / 2, 0.0f) - b_box->XMV_position;
	b_box->v_corner.push_back(c_corner);
	c_corner.XMV_position = XMVectorSet(XMF3_tmp.x + b_box->f_lengthX / 2, XMF3_tmp.y - b_box->f_lengthY / 2, XMF3_tmp.z - b_box->f_lengthZ / 2, 0.0f) - b_box->XMV_position;
	b_box->v_corner.push_back(c_corner);
	c_corner.XMV_position = XMVectorSet(XMF3_tmp.x + b_box->f_lengthX / 2, XMF3_tmp.y - b_box->f_lengthY / 2, XMF3_tmp.z + b_box->f_lengthZ / 2, 0.0f) - b_box->XMV_position;
	b_box->v_corner.push_back(c_corner);
	c_corner.XMV_position = XMVectorSet(XMF3_tmp.x + b_box->f_lengthX / 2, XMF3_tmp.y + b_box->f_lengthY / 2, XMF3_tmp.z - b_box->f_lengthZ / 2, 0.0f) - b_box->XMV_position;
	b_box->v_corner.push_back(c_corner);
	c_corner.XMV_position = XMVectorSet(XMF3_tmp.x + b_box->f_lengthX / 2, XMF3_tmp.y + b_box->f_lengthY / 2, XMF3_tmp.z + b_box->f_lengthZ / 2, 0.0f) - b_box->XMV_position;
	b_box->v_corner.push_back(c_corner);
}

void InitRBS()
{
	v_box.clear();
	XMFLOAT3X3 XMF3X3_inertiaTensor;
	XMMATRIX XMM_rotation = XMMatrixRotationZ(0.0f);
	XMFLOAT4 XMF4_tmp;
	XMMATRIX XMM_C = XMM_zero;
	XMFLOAT4X4 XMF4X4_C;

	//box1
	Box b_box(XMVectorSet(0.2f, 1.0f, 0.0f, 0.0f), 1.0f, 0.2f, 0.2f, 0.2f);
	b_box.XMV_orientation = XMQuaternionNormalize(XMQuaternionRotationMatrix(XMM_rotation));

	b_box.XMV_linearVelocity = XMV_zero;
	b_box.XMV_angularMomentum = XMV_zero;
	b_box.XMV_forceAccumulator = XMV_zero;
	b_box.XMV_torqueAccumulator = XMV_zero;

	CalculateCorners(&b_box);
	for (int i = 0; i < b_box.v_corner.size(); i++)
	{
		XMM_C += b_box.f_mass / b_box.v_corner.size() * VectorMulTranspose(b_box.v_corner[i].XMV_position);
	}

	XMStoreFloat4x4(&XMF4X4_C, XMM_C);
	float f_c = XMF4X4_C._11 + XMF4X4_C._22 + XMF4X4_C._33;
	b_box.XMM_inertiaTensor = XMMatrixIdentity() * f_c - XMM_C;

	b_box.XMM_inertiaTensor = XMMatrixRotationQuaternion(b_box.XMV_orientation) * XMMatrixInverse(NULL, b_box.XMM_inertiaTensor) * XMMatrixTranspose(XMMatrixRotationQuaternion(b_box.XMV_orientation));
	XMStoreFloat4x4(&XMF4X4_C, b_box.XMM_inertiaTensor);
	XMF4X4_C._44 = 0.0f;
	b_box.XMM_inertiaTensor = XMLoadFloat4x4(&XMF4X4_C);

	for (int i = 0; i < b_box.v_corner.size(); i++)
	{
		b_box.XMV_angularMomentum += XMVector3Cross(b_box.v_corner[i].XMV_position, (b_box.f_mass / b_box.v_corner.size()) * /*XMV_zero);//*/ b_box.v_corner[i].XMV_velocity);
	}

	b_box.XMV_angularVelocity = XMVector4Transform(b_box.XMV_angularMomentum, b_box.XMM_inertiaTensor);
	v_box.push_back(b_box);

	//box2
	XMM_rotation = XMMatrixRotationZ(0.0f);
	b_box = Box(XMVectorSet(0.0f, -0.5f + 0.05f, 0.0f, 0.0f), 1.0f, 0.3f, 0.2f, 0.1f);
	b_box.XMV_orientation = XMQuaternionNormalize(XMQuaternionRotationMatrix(XMM_rotation));

	b_box.XMV_linearVelocity = XMV_zero;
	b_box.XMV_angularMomentum = XMV_zero;
	b_box.XMV_forceAccumulator = XMV_zero;
	b_box.XMV_torqueAccumulator = XMV_zero;

	CalculateCorners(&b_box);
	for (int i = 0; i < b_box.v_corner.size(); i++)
	{
		XMM_C += b_box.f_mass / b_box.v_corner.size() * VectorMulTranspose(b_box.v_corner[i].XMV_position);
	}

	XMStoreFloat4x4(&XMF4X4_C, XMM_C);
	f_c = XMF4X4_C._11 + XMF4X4_C._22 + XMF4X4_C._33;
	b_box.XMM_inertiaTensor = XMMatrixIdentity() * f_c - XMM_C;

	b_box.XMM_inertiaTensor = XMMatrixRotationQuaternion(b_box.XMV_orientation) * XMMatrixInverse(NULL, b_box.XMM_inertiaTensor) * XMMatrixTranspose(XMMatrixRotationQuaternion(b_box.XMV_orientation));
	XMStoreFloat4x4(&XMF4X4_C, b_box.XMM_inertiaTensor);
	XMF4X4_C._44 = 0.0f;
	b_box.XMM_inertiaTensor = XMLoadFloat4x4(&XMF4X4_C);

	for (int i = 0; i < b_box.v_corner.size(); i++)
	{
		b_box.XMV_angularMomentum += XMVector3Cross(b_box.v_corner[i].XMV_position, (b_box.f_mass / b_box.v_corner.size()) * b_box.v_corner[i].XMV_velocity);
	}

	b_box.XMV_angularVelocity = XMVector4Transform(b_box.XMV_angularMomentum, b_box.XMM_inertiaTensor);
	v_box.push_back(b_box);
}

float CalculateImpulse(Box* b_boxA, Box* b_boxB, CollisionInfo col_info)
{
	XMVECTOR XMV_velPointA = b_boxA->XMV_linearVelocity + XMVector3Cross(b_boxA->XMV_angularVelocity, col_info.collisionPointWorld);
	//PrintVector(b_boxA->XMV_linearVelocity, "b_boxA->XMV_linearVelocity");
	//PrintVector(b_boxA->XMV_angularVelocity, "b_boxA->XMV_angularVelocity");
	//PrintVector(col_info.collisionPointWorld, "col_info.collisionPointWorld");

	//PrintVector(XMV_velPointA, "VelA");
	XMVECTOR XMV_velPointB = b_boxB->XMV_linearVelocity + XMVector3Cross(b_boxB->XMV_angularVelocity, col_info.collisionPointWorld);
	//PrintVector(XMV_velPointB, "VelB");
	XMVECTOR XMV_velRel = XMV_velPointA - XMV_velPointB;
	//PrintVector(XMV_velRel, "RelVel");

	XMVECTOR XMV_inertiaInfoA = XMVector3Cross(XMVector4Transform(XMVector3Cross(b_boxA->XMV_position, col_info.normalWorld), b_boxA->XMM_inertiaTensor), b_boxA->XMV_position);
	XMVECTOR XMV_inertiaInfoB = XMVector3Cross(XMVector4Transform(XMVector3Cross(b_boxB->XMV_position, col_info.normalWorld), b_boxB->XMM_inertiaTensor), b_boxB->XMV_position);

	XMVECTOR XMV_dot1 = XMVector3Dot(XMV_velRel, col_info.normalWorld);

	XMVECTOR XMV_dot2 = XMVector3Dot((XMV_inertiaInfoA + XMV_inertiaInfoB), col_info.normalWorld);
	//PrintVector(XMV_dot1, "Dot1");
	//PrintVector(XMV_dot2, "Dot2");


	float f_impulse = -((1 + g_fImpulseConstant) * XMVectorGetX(XMV_dot1) / (1 / b_boxA->f_mass + 1 / b_boxB->f_mass + XMVectorGetX(XMV_dot2)));
	//std::cout << "Impulse: " << f_impulse << "\n";
	return f_impulse;
}

void CollisionDetectionRigidbody()
{
	for (int i = 0; i < v_box.size(); i++)
	{
		XMFLOAT3 tmp;
		XMStoreFloat3(&tmp, v_box[i].XMV_position);
		if (tmp.y < -0.5f + v_box[i].f_lengthY / 2.0f)
		{
			XMStoreFloat3(&tmp, v_box[i].XMV_linearVelocity);
			tmp.y = 0.0f;// -f_gravity * v_box[i].f_mass / v_box[i].v_corner.size();
			v_box[i].XMV_linearVelocity = XMLoadFloat3(&tmp);
		}
	}
	XMFLOAT4 XMF4_tmp1;
	XMFLOAT4 XMF4_tmp2;
	CollisionInfo col_info;
	for (int i = 0; i < v_box.size(); i++)
	{
		for (int k = 0; k < v_box.size(); k++)
		{
			if (i != k)
			{
				XMStoreFloat4(&XMF4_tmp1, v_box[i].XMV_position);
				XMStoreFloat4(&XMF4_tmp2, v_box[k].XMV_position);
				col_info = checkCollision(XMMatrixTranslation(XMF4_tmp1.x, XMF4_tmp1.y, XMF4_tmp1.z), XMMatrixTranslation(XMF4_tmp2.x, XMF4_tmp2.y, XMF4_tmp2.z), v_box[i].f_lengthX, v_box[i].f_lengthY, v_box[i].f_lengthZ, v_box[k].f_lengthX, v_box[k].f_lengthY, v_box[k].f_lengthZ);
				if (col_info.isValid)
				{
					XMFLOAT3 XMF3_hitPoint;
					XMStoreFloat3(&XMF3_hitPoint, col_info.collisionPointWorld);
					//std::cout << "hit! in " << XMF3_hitPoint.x << ":" << XMF3_hitPoint.y << ":" << XMF3_hitPoint.z << "\n";

					float f_impulse = CalculateImpulse(&v_box[i], &v_box[k], col_info);
					v_box[i].XMV_linearVelocity = v_box[i].XMV_linearVelocity + f_impulse * col_info.normalWorld / v_box[i].f_mass;
					v_box[k].XMV_linearVelocity = v_box[k].XMV_linearVelocity - f_impulse * col_info.normalWorld / v_box[k].f_mass;

					v_box[i].XMV_angularMomentum = v_box[i].XMV_angularMomentum + XMVector3Cross(v_box[i].XMV_position, f_impulse * col_info.normalWorld);
					v_box[k].XMV_angularMomentum = v_box[k].XMV_angularMomentum - XMVector3Cross(v_box[k].XMV_position, f_impulse * col_info.normalWorld);
				}
			}
		}
	}
}

void ApplyGravityToRigidbodys()
{
	XMVECTOR XMV_gravity = XMVectorSet(0.0f, -9.81f, 0.0f, 0.0f);
	for (int i = 0; i < v_box.size(); i++)
	{
		v_box[i].XMV_position += XMV_gravity * g_fTimeStepSize * 0.01f;
		for (int k = 0; k < v_box[i].v_corner.size(); k++)
		{
			v_box[i].v_corner[k].XMV_velocity = XMV_gravity * g_fTimeStepSize;
		}
		//XMFLOAT4X4 tmp4X4;
		//XMFLOAT3 tmp3;
		//XMStoreFloat3(&tmp3, v_box[i].XMV_angularMomentum);
		//XMStoreFloat4x4(&tmp4X4, v_box[i].XMM_inertiaTensor);
		//std::cout << "angularMomentum: " << tmp3.x << " " << tmp3.y << " " << tmp3.z << "\n";
		//std::cout << "inertiaTensor: " << tmp4X4._11 << " " << tmp4X4._22 << " " << tmp4X4._33 << " " << tmp4X4._44 << "\n";
		//XMStoreFloat3(&tmp3, v_box[i].XMV_angularVelocity);
		//std::cout << "angularVelocity: " << tmp3.x << " " << tmp3.y << " " << tmp3.z << "\n";
	}
}

void ApplyPhysikRBS()
{
	for (int i = 0; i < v_box.size(); i++)
	{
		v_box[i].XMV_forceAccumulator = XMV_zero;
		v_box[i].XMV_torqueAccumulator = XMV_zero;
		for (int k = 0; k < v_box[i].v_corner.size(); k++)
		{
			v_box[i].XMV_forceAccumulator += v_box[i].v_corner[k].XMV_force;
			v_box[i].XMV_torqueAccumulator += XMVector3Cross(v_box[i].v_corner[k].XMV_position, v_box[i].v_corner[k].XMV_force);
		}
		//PrintVector(v_box[0].XMV_torqueAccumulator, "XMV_torqueAccumulator");

		v_box[i].XMV_position += v_box[i].XMV_linearVelocity * g_fTimeStepSize;
		v_box[i].XMV_linearVelocity += (v_box[i].XMV_forceAccumulator * g_fTimeStepSize) / v_box[i].f_mass;
		//PrintVector(v_box[0].XMV_linearVelocity, "XMV_linearVelocity");
		//PrintVector((v_box[i].XMV_forceAccumulator * g_fTimeStepSize) / v_box[i].f_mass, "Berechnung");

		v_box[i].XMV_orientation = g_fTimeStepSize / 2 * XMQuaternionMultiply(VectorToQuaternion(v_box[i].XMV_angularVelocity), v_box[i].XMV_orientation);
		v_box[i].XMV_orientation = XMQuaternionNormalize(v_box[i].XMV_orientation);
		v_box[i].XMV_angularMomentum += v_box[i].XMV_torqueAccumulator * g_fTimeStepSize;

		v_box[i].XMM_inertiaTensor = XMMatrixRotationQuaternion(v_box[i].XMV_orientation) * v_box[i].XMM_inertiaTensor * XMMatrixTranspose(XMMatrixRotationQuaternion(v_box[i].XMV_orientation));

		v_box[i].XMV_angularVelocity += XMVector4Transform(v_box[i].XMV_angularMomentum, v_box[i].XMM_inertiaTensor);

		//PrintVector(v_box[0].XMV_angularVelocity, "XMV_angularVelocity");
		for (int k = 0; k < v_box[i].v_corner.size(); k++)
		{
			v_box[i].v_corner[k].XMV_position = v_box[i].XMV_position + XMVector3Transform(v_box[i].v_corner[k].XMV_position, XMMatrixRotationQuaternion(v_box[i].XMV_orientation));
			v_box[i].v_corner[k].XMV_velocity = v_box[i].XMV_linearVelocity + XMVector3Cross(v_box[i].XMV_angularVelocity, v_box[i].v_corner[k].XMV_position);
		}
	}
}



//////////////MASSCOLLISIONDETECTION/////////////////
void InitBalls()
{
	v_ball.clear();
	Ball ba_ball;
	float f_distance = g_fSphereSize * 2;
	int i_count = 0;
	for (int i = 0; i < g_iNumSpheres; i++)
	{
		for (int k = 0; k < g_iNumSpheres; k++)
		{
			for (int l = 0; l < g_iNumSpheres; l++)
			{
				ba_ball.XMV_position = XMVectorSet(-0.25f + i * f_distance, -0.4f + k * f_distance, -0.25f + l * f_distance, 0.0f);
				ba_ball.XMV_velocity = XMV_zero;
				ba_ball.id = i_count;
				v_ball.push_back(ba_ball);
				i_count++;
			}
		}
	}
}

//Restraining Balls to th drawn cube(0.5, 0.5, 0.5)
void RestrainingPosition()
{
	for (int i = 0; i < v_ball.size(); i++)
	{
		//ToDo: change in velocoty instead of position
		if (XMVectorGetX(v_ball[i].XMV_position) < -0.5f)
			v_ball[i].XMV_velocity = XMVectorScale(XMVectorSetX(v_ball[i].XMV_velocity, -XMVectorGetX(v_ball[i].XMV_velocity)), g_fDamping);
		if (XMVectorGetX(v_ball[i].XMV_position) > 0.5f)
			v_ball[i].XMV_velocity = XMVectorScale(XMVectorSetX(v_ball[i].XMV_velocity, -XMVectorGetX(v_ball[i].XMV_velocity)), g_fDamping);
		if (XMVectorGetY(v_ball[i].XMV_position) < -0.5f)
			v_ball[i].XMV_velocity = XMVectorScale(XMVectorSetY(v_ball[i].XMV_velocity, -XMVectorGetY(v_ball[i].XMV_velocity)), g_fDamping);
		if (XMVectorGetY(v_ball[i].XMV_position) > 0.5f)
			v_ball[i].XMV_velocity = XMVectorScale(XMVectorSetY(v_ball[i].XMV_velocity, -XMVectorGetY(v_ball[i].XMV_velocity)), g_fDamping);
		if (XMVectorGetZ(v_ball[i].XMV_position) < -0.5f)
			v_ball[i].XMV_velocity = XMVectorScale(XMVectorSetZ(v_ball[i].XMV_velocity, -XMVectorGetZ(v_ball[i].XMV_velocity)), g_fDamping);
		if (XMVectorGetZ(v_ball[i].XMV_position) > 0.5f)
			v_ball[i].XMV_velocity = XMVectorScale(XMVectorSetZ(v_ball[i].XMV_velocity, -XMVectorGetZ(v_ball[i].XMV_velocity)), g_fDamping);
	}
}

//Apllys gravity to balls
void ApplyGravityToRigidbodysToBalls()
{
	for (int i = 0; i < v_ball.size(); i++)
	{
		v_ball[i].XMV_velocity = XMVectorAdd(v_ball[i].XMV_velocity, XMVectorScale(XMVectorScale(XMVectorScale(XMVectorSet(0.0f, f_gravity, 0.0f, 0.0f), 1 / g_fMass), g_fTimeStepSize), g_fDamping));
			//XMVectorSet(0.0f, f_gravity, 0.0f, 0.0f) / g_fMass * g_fTimeStepSize * g_fDamping;
			//XMVectorScale(XMVectorScale(XMVectorScale(XMVectorAdd(v_ball[i].XMV_velocity, XMVectorSet(0.0f, f_gravity, 0.0f, 0.0f)), 1 / g_fMass), g_fTimeStepSize), g_fDamping);
	}
}

//Computes Impuls of BallCollision
void BallCollisionImpuls(Ball* ba_sphereA, Ball* ba_sphereB)
{
	XMVECTOR XMV_normal = XMVectorSubtract(ba_sphereA->XMV_position, ba_sphereB->XMV_position);
	float f_impuls = g_fcollisionScalar * (1 - XMVectorGetX(XMVector3Length(XMV_normal)) / (g_fSphereSize * 2));
	//std::cout << f_impuls << "\n";
	ba_sphereA->XMV_velocity = XMVectorAdd(ba_sphereA->XMV_velocity, XMVectorScale(XMV_normal,1 / g_fMass * g_fDamping * f_impuls));
	ba_sphereB->XMV_velocity = XMVectorAdd(ba_sphereB->XMV_velocity, XMVectorScale(XMV_normal, 1 / g_fMass * g_fDamping * -f_impuls));
}

void NaiveCollisionDetection()
{
	for (int i = 0; i < v_ball.size(); i++)
	{
		// k = i + 1
		for (int k = 0; k < v_ball.size(); k++)
		{
			if (XMVectorGetX(XMVector3Length(XMVectorSubtract(v_ball[i].XMV_position, v_ball[k].XMV_position))) < g_fSphereSize * 2)
			{
				BallCollisionImpuls(&v_ball[i], &v_ball[k]);
			}
		}
		RestrainingPosition();
	}
}

void UpdateBallPosition()
{
	for (int i = 0; i < v_ball.size(); i++)
	{
		v_ball[i].XMV_position = XMVectorAdd(v_ball[i].XMV_position, v_ball[i].XMV_velocity);
	}
}



// Draw the edges of the bounding box [-0.5;0.5]� rotated with the cameras model tranformation.
// (Drawn as line primitives using a DirectXTK primitive batch)
void DrawBoundingBox(ID3D11DeviceContext* pd3dImmediateContext)
{
	// Setup position/color effect
	g_pEffectPositionColor->SetWorld(g_camera.GetWorldMatrix());

	g_pEffectPositionColor->Apply(pd3dImmediateContext);
	pd3dImmediateContext->IASetInputLayout(g_pInputLayoutPositionColor);

	// Draw
	g_pPrimitiveBatchPositionColor->Begin();

	// Lines in x direction (red color)
	for (int i = 0; i < 4; i++)
	{
		g_pPrimitiveBatchPositionColor->DrawLine(
			VertexPositionColor(XMVectorSet(-0.5f, (float)(i % 2) - 0.5f, (float)(i / 2) - 0.5f, 1), Colors::Red),
			VertexPositionColor(XMVectorSet(0.5f, (float)(i % 2) - 0.5f, (float)(i / 2) - 0.5f, 1), Colors::Red)
			);
	}

	// Lines in y direction
	for (int i = 0; i < 4; i++)
	{
		g_pPrimitiveBatchPositionColor->DrawLine(
			VertexPositionColor(XMVectorSet((float)(i % 2) - 0.5f, -0.5f, (float)(i / 2) - 0.5f, 1), Colors::Green),
			VertexPositionColor(XMVectorSet((float)(i % 2) - 0.5f, 0.5f, (float)(i / 2) - 0.5f, 1), Colors::Green)
			);
	}

	// Lines in z direction
	for (int i = 0; i < 4; i++)
	{
		g_pPrimitiveBatchPositionColor->DrawLine(
			VertexPositionColor(XMVectorSet((float)(i % 2) - 0.5f, (float)(i / 2) - 0.5f, -0.5f, 1), Colors::Blue),
			VertexPositionColor(XMVectorSet((float)(i % 2) - 0.5f, (float)(i / 2) - 0.5f, 0.5f, 1), Colors::Blue)
			);
	}

	g_pPrimitiveBatchPositionColor->End();
}

// Draw a large, square plane at y=-1 with a checkerboard pattern
// (Drawn as multiple quads, i.e. triangle strips, using a DirectXTK primitive batch)
void DrawFloor(ID3D11DeviceContext* pd3dImmediateContext)
{
	// Setup position/normal/color effect
	g_pEffectPositionNormalColor->SetWorld(XMMatrixIdentity());
	g_pEffectPositionNormalColor->SetEmissiveColor(Colors::Black);
	g_pEffectPositionNormalColor->SetDiffuseColor(0.8f * Colors::White);
	g_pEffectPositionNormalColor->SetSpecularColor(0.4f * Colors::White);
	g_pEffectPositionNormalColor->SetSpecularPower(1000);

	g_pEffectPositionNormalColor->Apply(pd3dImmediateContext);
	pd3dImmediateContext->IASetInputLayout(g_pInputLayoutPositionNormalColor);

	// Draw 4*n*n quads spanning x = [-n;n], y = -1, z = [-n;n]
	const float n = 4;
	XMVECTOR normal = XMVectorSet(0, 1, 0, 0);
	XMVECTOR planecenter = XMVectorSet(0, -2, 0, 0);

	g_pPrimitiveBatchPositionNormalColor->Begin();
	for (float z = -n; z < n; z++)
	{
		for (float x = -n; x < n; x++)
		{
			// Quad vertex positions
			XMVECTOR pos[] = { XMVectorSet(x, -1, z + 1, 0),
				XMVectorSet(x + 1, -1, z + 1, 0),
				XMVectorSet(x + 1, -1, z, 0),
				XMVectorSet(x, -1, z, 0) };

			// Color checkerboard pattern (white & gray)
			XMVECTOR color = ((int(z + x) % 2) == 0) ? XMVectorSet(1, 1, 1, 1) : XMVectorSet(0.6f, 0.6f, 0.6f, 1);

			// Color attenuation based on distance to plane center
			float attenuation[] = {
				1.0f - XMVectorGetX(XMVector3Length(pos[0] - planecenter)) / n,
				1.0f - XMVectorGetX(XMVector3Length(pos[1] - planecenter)) / n,
				1.0f - XMVectorGetX(XMVector3Length(pos[2] - planecenter)) / n,
				1.0f - XMVectorGetX(XMVector3Length(pos[3] - planecenter)) / n };

			g_pPrimitiveBatchPositionNormalColor->DrawQuad(
				VertexPositionNormalColor(pos[0], normal, attenuation[0] * color),
				VertexPositionNormalColor(pos[1], normal, attenuation[1] * color),
				VertexPositionNormalColor(pos[2], normal, attenuation[2] * color),
				VertexPositionNormalColor(pos[3], normal, attenuation[3] * color)
				);
		}
	}
	g_pPrimitiveBatchPositionNormalColor->End();
}

// Draw several objects randomly positioned in [-0.5f;0.5]�  using DirectXTK geometric primitives.
void DrawSomeRandomObjects(ID3D11DeviceContext* pd3dImmediateContext)
{
	// Setup position/normal effect (constant variables)
	g_pEffectPositionNormal->SetEmissiveColor(Colors::Black);
	g_pEffectPositionNormal->SetSpecularColor(0.4f * Colors::White);
	g_pEffectPositionNormal->SetSpecularPower(100);

	std::mt19937 eng;
	std::uniform_real_distribution<float> randCol(0.0f, 1.0f);
	std::uniform_real_distribution<float> randPos(-0.5f, 0.5f);

	for (int i = 0; i < g_iNumSpheres; i++)
	{
		// Setup position/normal effect (per object variables)
		g_pEffectPositionNormal->SetDiffuseColor(0.6f * XMColorHSVToRGB(XMVectorSet(randCol(eng), 1, 1, 0)));
		XMMATRIX scale = XMMatrixScaling(g_fSphereSize, g_fSphereSize, g_fSphereSize);
		XMMATRIX trans = XMMatrixTranslation(randPos(eng), randPos(eng), randPos(eng));
		g_pEffectPositionNormal->SetWorld(scale * trans * g_camera.GetWorldMatrix());

		// Draw
		// NOTE: The following generates one draw call per object, so performance will be bad for n>>1000 or so
		g_pSphere->Draw(g_pEffectPositionNormal, g_pInputLayoutPositionNormal);
	}
}

// Draw a teapot at the position g_vfMovableObjectPos.
void DrawMovableTeapot(ID3D11DeviceContext* pd3dImmediateContext)
{
	// Setup position/normal effect (constant variables)
	g_pEffectPositionNormal->SetEmissiveColor(Colors::Black);
	g_pEffectPositionNormal->SetDiffuseColor(0.6f * Colors::Cornsilk);
	g_pEffectPositionNormal->SetSpecularColor(0.4f * Colors::White);
	g_pEffectPositionNormal->SetSpecularPower(100);

	XMMATRIX scale = XMMatrixScaling(0.5f, 0.5f, 0.5f);
	XMMATRIX trans = XMMatrixTranslation(g_vfMovableObjectPos.x, g_vfMovableObjectPos.y, g_vfMovableObjectPos.z);
	g_pEffectPositionNormal->SetWorld(scale * trans);

	// Draw
	g_pTeapot->Draw(g_pEffectPositionNormal, g_pInputLayoutPositionNormal);
}

// Draw a simple triangle using custom shaders (g_pEffect)
void DrawTriangle(ID3D11DeviceContext* pd3dImmediateContext)
{
	XMMATRIX world = g_camera.GetWorldMatrix();
	XMMATRIX view = g_camera.GetViewMatrix();
	XMMATRIX proj = g_camera.GetProjMatrix();
	XMFLOAT4X4 mViewProj;
	XMStoreFloat4x4(&mViewProj, world * view * proj);
	g_pEffect->GetVariableByName("g_worldViewProj")->AsMatrix()->SetMatrix((float*)mViewProj.m);
	g_pEffect->GetTechniqueByIndex(0)->GetPassByIndex(0)->Apply(0, pd3dImmediateContext);

	pd3dImmediateContext->IASetVertexBuffers(0, 0, nullptr, nullptr, nullptr);
	pd3dImmediateContext->IASetIndexBuffer(nullptr, DXGI_FORMAT_R16_UINT, 0);
	pd3dImmediateContext->IASetInputLayout(nullptr);
	pd3dImmediateContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
	pd3dImmediateContext->Draw(3, 0);
}

void DrawMassSpringSystem(ID3D11DeviceContext* pd3dImmediateContext)
{
	g_pEffectPositionNormal->SetEmissiveColor(Colors::Black);
	g_pEffectPositionNormal->SetSpecularColor(0.4f * Colors::White);
	g_pEffectPositionNormal->SetSpecularPower(100);

	if (b_start)
	{
		InitPoints();
		InitSprings();
		b_start = false;
		TwDeleteBar(g_pTweakBar);
		g_pTweakBar = nullptr;
		TwTerminate();
		InitTweakBar(pd3dDeviceTweakbar);
	}

	for (int i = 0; i < v_point.size(); i++)
	{
		//set position in worldspace
		XMMATRIX XMM_scale = XMMatrixScaling(g_fSphereSize, g_fSphereSize, g_fSphereSize);
		XMMATRIX XMM_trans = XMMatrixTranslationFromVector(v_point[i].XMV_position);
		XMMATRIX XMM_sphereWorldMatrix = XMM_scale * XMM_trans * g_camera.GetWorldMatrix();

		g_pEffectPositionNormal->SetWorld(XMM_sphereWorldMatrix);

		if (v_point[i].b_Dummy)
			g_pEffectPositionNormal->SetDiffuseColor(0.6f * XMColorHSVToRGB(XMVectorSet(0.5, 1, 1, 0)));
		else
			g_pEffectPositionNormal->SetDiffuseColor(0.6f * XMColorHSVToRGB(XMVectorSet(1, 1, 1, 0)));
		g_pSphere->Draw(g_pEffectPositionNormal, g_pInputLayoutPositionNormal);
	}

	g_pEffectPositionColor->SetWorld(g_camera.GetWorldMatrix());

	g_pEffectPositionColor->Apply(pd3dImmediateContext);
	pd3dImmediateContext->IASetInputLayout(g_pInputLayoutPositionColor);
	g_pPrimitiveBatchPositionColor->Begin();

	for (int i = 0; i < v_spring.size(); i++){
		g_pPrimitiveBatchPositionColor->DrawLine(
			VertexPositionColor(GetPointOf(v_spring[i].i_point1).XMV_position, Colors::Green),
			VertexPositionColor(GetPointOf(v_spring[i].i_point2).XMV_position, Colors::Green));
	}
	g_pPrimitiveBatchPositionColor->End();
}

void DrawRigidBody(ID3D11DeviceContext* pd3dImmediateContext)
{
	//ToDo
	// Setup position/normal effect (constant variables)
	g_pEffectPositionNormal->SetEmissiveColor(Colors::Black);
	g_pEffectPositionNormal->SetDiffuseColor(0.6f * Colors::Cornsilk);
	g_pEffectPositionNormal->SetSpecularColor(0.4f * Colors::White);
	g_pEffectPositionNormal->SetSpecularPower(100);

	if (b_start)
	{
		InitRBS();
		b_start = false;
		TwDeleteBar(g_pTweakBar);
		g_pTweakBar = nullptr;
		TwTerminate();
		InitTweakBar(pd3dDeviceTweakbar);
	}

	for (int i = 0; i < v_box.size(); i++)
	{
		XMFLOAT3 tmp;
		XMStoreFloat3(&tmp, v_box[i].XMV_position);
		XMMATRIX scale = XMMatrixScaling(v_box[i].f_lengthX * 2, v_box[i].f_lengthY * 2, v_box[i].f_lengthZ * 2);
		XMMATRIX trans = XMMatrixTranslation(tmp.x, tmp.y, tmp.z);
		XMMATRIX rotate = XMMatrixRotationQuaternion(v_box[i].XMV_orientation);
		g_pEffectPositionNormal->SetWorld(scale * rotate *trans);
		g_pCube->Draw(g_pEffectPositionNormal, g_pInputLayoutPositionNormal);
	}

	//// Setup position/color effect
	//g_pEffectPositionColor->SetWorld(g_camera.GetWorldMatrix());
	//
	//g_pEffectPositionColor->Apply(pd3dImmediateContext);
	//pd3dImmediateContext->IASetInputLayout(g_pInputLayoutPositionColor);
	//
	//// Draw
	//g_pPrimitiveBatchPositionColor->Begin();
	//
	//for (int i = 0; i < v_box.size(); i++)
	//{
	//	XMFLOAT3 tmp;
	//	XMStoreFloat3(&tmp, v_box[i].XMV_position);
	//
	//	// Lines in x direction (red color)
	//	g_pPrimitiveBatchPositionColor->DrawLine(
	//		VertexPositionColor(XMVectorSet(tmp.x - v_box[i].f_lengthX / 2.0f,tmp.y - v_box[i].f_lengthY / 2.0f, tmp.z + v_box[i].f_lengthZ / 2.0f, 1), Colors::Red),
	//		VertexPositionColor(XMVectorSet(tmp.x + v_box[i].f_lengthX / 2.0f,tmp.y - v_box[i].f_lengthY / 2.0f, tmp.z + v_box[i].f_lengthZ / 2.0f, 1), Colors::Red)
	//		);
	//	g_pPrimitiveBatchPositionColor->DrawLine(
	//		VertexPositionColor(XMVectorSet(tmp.x - v_box[i].f_lengthX / 2.0f, tmp.y + v_box[i].f_lengthY / 2.0f, tmp.z + v_box[i].f_lengthZ / 2.0f, 1), Colors::Red),
	//		VertexPositionColor(XMVectorSet(tmp.x + v_box[i].f_lengthX / 2.0f, tmp.y + v_box[i].f_lengthY / 2.0f, tmp.z + v_box[i].f_lengthZ / 2.0f, 1), Colors::Red)
	//		);
	//	g_pPrimitiveBatchPositionColor->DrawLine(
	//		VertexPositionColor(XMVectorSet(tmp.x - v_box[i].f_lengthX / 2.0f, tmp.y - v_box[i].f_lengthY / 2.0f, tmp.z - v_box[i].f_lengthZ / 2.0f, 1), Colors::Red),
	//		VertexPositionColor(XMVectorSet(tmp.x + v_box[i].f_lengthX / 2.0f, tmp.y - v_box[i].f_lengthY / 2.0f, tmp.z - v_box[i].f_lengthZ / 2.0f, 1), Colors::Red)
	//		);
	//	g_pPrimitiveBatchPositionColor->DrawLine(
	//		VertexPositionColor(XMVectorSet(tmp.x - v_box[i].f_lengthX / 2.0f, tmp.y + v_box[i].f_lengthY / 2.0f, tmp.z - v_box[i].f_lengthZ / 2.0f, 1), Colors::Red),
	//		VertexPositionColor(XMVectorSet(tmp.x + v_box[i].f_lengthX / 2.0f, tmp.y + v_box[i].f_lengthY / 2.0f, tmp.z - v_box[i].f_lengthZ / 2.0f, 1), Colors::Red)
	//		);
	//
	//	// Lines in y direction
	//	g_pPrimitiveBatchPositionColor->DrawLine(
	//		VertexPositionColor(XMVectorSet(tmp.x + v_box[i].f_lengthX / 2.0f, tmp.y - v_box[i].f_lengthY / 2.0f, tmp.z + v_box[i].f_lengthZ / 2.0f, 1), Colors::Red),
	//		VertexPositionColor(XMVectorSet(tmp.x + v_box[i].f_lengthX / 2.0f, tmp.y + v_box[i].f_lengthY / 2.0f, tmp.z + v_box[i].f_lengthZ / 2.0f, 1), Colors::Red)
	//		);
	//	g_pPrimitiveBatchPositionColor->DrawLine(
	//		VertexPositionColor(XMVectorSet(tmp.x + v_box[i].f_lengthX / 2.0f, tmp.y - v_box[i].f_lengthY / 2.0f, tmp.z - v_box[i].f_lengthZ / 2.0f, 1), Colors::Red),
	//		VertexPositionColor(XMVectorSet(tmp.x + v_box[i].f_lengthX / 2.0f, tmp.y + v_box[i].f_lengthY / 2.0f, tmp.z - v_box[i].f_lengthZ / 2.0f, 1), Colors::Red)
	//		);
	//	g_pPrimitiveBatchPositionColor->DrawLine(
	//		VertexPositionColor(XMVectorSet(tmp.x - v_box[i].f_lengthX / 2.0f, tmp.y - v_box[i].f_lengthY / 2.0f, tmp.z + v_box[i].f_lengthZ / 2.0f, 1), Colors::Red),
	//		VertexPositionColor(XMVectorSet(tmp.x - v_box[i].f_lengthX / 2.0f, tmp.y + v_box[i].f_lengthY / 2.0f, tmp.z + v_box[i].f_lengthZ / 2.0f, 1), Colors::Red)
	//		);
	//	g_pPrimitiveBatchPositionColor->DrawLine(
	//		VertexPositionColor(XMVectorSet(tmp.x - v_box[i].f_lengthX / 2.0f, tmp.y - v_box[i].f_lengthY / 2.0f, tmp.z - v_box[i].f_lengthZ / 2.0f, 1), Colors::Red),
	//		VertexPositionColor(XMVectorSet(tmp.x - v_box[i].f_lengthX / 2.0f, tmp.y + v_box[i].f_lengthY / 2.0f, tmp.z - v_box[i].f_lengthZ / 2.0f, 1), Colors::Red)
	//		);
	//	// Lines in z direction
	//	g_pPrimitiveBatchPositionColor->DrawLine(
	//		VertexPositionColor(XMVectorSet(tmp.x + v_box[i].f_lengthX / 2.0f, tmp.y + v_box[i].f_lengthY / 2.0f, tmp.z - v_box[i].f_lengthZ / 2.0f, 1), Colors::Red),
	//		VertexPositionColor(XMVectorSet(tmp.x + v_box[i].f_lengthX / 2.0f, tmp.y + v_box[i].f_lengthY / 2.0f, tmp.z + v_box[i].f_lengthZ / 2.0f, 1), Colors::Red)
	//		);
	//	g_pPrimitiveBatchPositionColor->DrawLine(
	//		VertexPositionColor(XMVectorSet(tmp.x + v_box[i].f_lengthX / 2.0f, tmp.y - v_box[i].f_lengthY / 2.0f, tmp.z - v_box[i].f_lengthZ / 2.0f, 1), Colors::Red),
	//		VertexPositionColor(XMVectorSet(tmp.x + v_box[i].f_lengthX / 2.0f, tmp.y - v_box[i].f_lengthY / 2.0f, tmp.z + v_box[i].f_lengthZ / 2.0f, 1), Colors::Red)
	//		);
	//	g_pPrimitiveBatchPositionColor->DrawLine(
	//		VertexPositionColor(XMVectorSet(tmp.x - v_box[i].f_lengthX / 2.0f, tmp.y + v_box[i].f_lengthY / 2.0f, tmp.z - v_box[i].f_lengthZ / 2.0f, 1), Colors::Red),
	//		VertexPositionColor(XMVectorSet(tmp.x - v_box[i].f_lengthX / 2.0f, tmp.y + v_box[i].f_lengthY / 2.0f, tmp.z + v_box[i].f_lengthZ / 2.0f, 1), Colors::Red)
	//		);
	//	g_pPrimitiveBatchPositionColor->DrawLine(
	//		VertexPositionColor(XMVectorSet(tmp.x - v_box[i].f_lengthX / 2.0f, tmp.y - v_box[i].f_lengthY / 2.0f, tmp.z - v_box[i].f_lengthZ / 2.0f, 1), Colors::Red),
	//		VertexPositionColor(XMVectorSet(tmp.x - v_box[i].f_lengthX / 2.0f, tmp.y - v_box[i].f_lengthY / 2.0f, tmp.z + v_box[i].f_lengthZ / 2.0f, 1), Colors::Red)
	//		);
	//}
	//	g_pPrimitiveBatchPositionColor->End();		
}

void DrawBiab(ID3D11DeviceContext* pd3dImmediateContext)
{
	g_pEffectPositionNormal->SetEmissiveColor(Colors::Black);
	g_pEffectPositionNormal->SetSpecularColor(0.4f * Colors::White);
	g_pEffectPositionNormal->SetSpecularPower(100);

	if (b_start)
	{
		InitBalls();
		b_start = false;
		TwDeleteBar(g_pTweakBar);
		g_pTweakBar = nullptr;
		TwTerminate();
		InitTweakBar(pd3dDeviceTweakbar);
	}

	for (int i = 0; i < v_ball.size(); i++)
	{
		//set position in worldspace
		XMMATRIX XMM_scale = XMMatrixScaling(g_fSphereSize, g_fSphereSize, g_fSphereSize);
		XMMATRIX XMM_trans = XMMatrixTranslationFromVector(v_ball[i].XMV_position);
		XMMATRIX XMM_sphereWorldMatrix = XMM_scale * XMM_trans * g_camera.GetWorldMatrix();

		g_pEffectPositionNormal->SetWorld(XMM_sphereWorldMatrix);
		g_pSphere->Draw(g_pEffectPositionNormal, g_pInputLayoutPositionNormal);
	}
}

// ============================================================
// DXUT Callbacks
// ============================================================


//--------------------------------------------------------------------------------------
// Reject any D3D11 devices that aren't acceptable by returning false
//--------------------------------------------------------------------------------------
bool CALLBACK IsD3D11DeviceAcceptable(const CD3D11EnumAdapterInfo *AdapterInfo, UINT Output, const CD3D11EnumDeviceInfo *DeviceInfo,
	DXGI_FORMAT BackBufferFormat, bool bWindowed, void* pUserContext)
{
	return true;
}


//--------------------------------------------------------------------------------------
// Called right before creating a device, allowing the app to modify the device settings as needed
//--------------------------------------------------------------------------------------
bool CALLBACK ModifyDeviceSettings(DXUTDeviceSettings* pDeviceSettings, void* pUserContext)
{
	return true;
}


//--------------------------------------------------------------------------------------
// Create any D3D11 resources that aren't dependent on the back buffer
//--------------------------------------------------------------------------------------
HRESULT CALLBACK OnD3D11CreateDevice(ID3D11Device* pd3dDevice, const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc, void* pUserContext)
{
	HRESULT hr;

	ID3D11DeviceContext* pd3dImmediateContext = DXUTGetD3D11DeviceContext();;

	std::wcout << L"Device: " << DXUTGetDeviceStats() << std::endl;

	// Load custom effect from "effect.fxo" (compiled "effect.fx")
	std::wstring effectPath = GetExePath() + L"effect.fxo";
	if (FAILED(hr = D3DX11CreateEffectFromFile(effectPath.c_str(), 0, pd3dDevice, &g_pEffect)))
	{
		std::wcout << L"Failed creating effect with error code " << int(hr) << std::endl;
		return hr;
	}

	pd3dDeviceTweakbar = pd3dDevice;

	// Init AntTweakBar GUI
	InitTweakBar(pd3dDevice);

	// Create DirectXTK geometric primitives for later usage
	g_pSphere = GeometricPrimitive::CreateGeoSphere(pd3dImmediateContext, 2.0f, 2, false);
	g_pTeapot = GeometricPrimitive::CreateTeapot(pd3dImmediateContext, 1.5f, 8, false);
	g_pCube = GeometricPrimitive::CreateCube(pd3dImmediateContext, 0.5f, false);

	// Create effect, input layout and primitive batch for position/color vertices (DirectXTK)
	{
		// Effect
		g_pEffectPositionColor = new BasicEffect(pd3dDevice);
		g_pEffectPositionColor->SetVertexColorEnabled(true); // triggers usage of position/color vertices

		// Input layout
		void const* shaderByteCode;
		size_t byteCodeLength;
		g_pEffectPositionColor->GetVertexShaderBytecode(&shaderByteCode, &byteCodeLength);

		pd3dDevice->CreateInputLayout(VertexPositionColor::InputElements,
			VertexPositionColor::InputElementCount,
			shaderByteCode, byteCodeLength,
			&g_pInputLayoutPositionColor);

		// Primitive batch
		g_pPrimitiveBatchPositionColor = new PrimitiveBatch<VertexPositionColor>(pd3dImmediateContext);
	}

	// Create effect, input layout and primitive batch for position/normal vertices (DirectXTK)
	{
		// Effect
		g_pEffectPositionNormal = new BasicEffect(pd3dDevice);
		g_pEffectPositionNormal->EnableDefaultLighting(); // triggers usage of position/normal vertices
		g_pEffectPositionNormal->SetPerPixelLighting(true);

		// Input layout
		void const* shaderByteCode;
		size_t byteCodeLength;
		g_pEffectPositionNormal->GetVertexShaderBytecode(&shaderByteCode, &byteCodeLength);

		pd3dDevice->CreateInputLayout(VertexPositionNormal::InputElements,
			VertexPositionNormal::InputElementCount,
			shaderByteCode, byteCodeLength,
			&g_pInputLayoutPositionNormal);

		// Primitive batch
		g_pPrimitiveBatchPositionNormal = new PrimitiveBatch<VertexPositionNormal>(pd3dImmediateContext);
	}

	// Create effect, input layout and primitive batch for position/normal/color vertices (DirectXTK)
	{
		// Effect
		g_pEffectPositionNormalColor = new BasicEffect(pd3dDevice);
		g_pEffectPositionNormalColor->SetPerPixelLighting(true);
		g_pEffectPositionNormalColor->EnableDefaultLighting();     // triggers usage of position/normal/color vertices
		g_pEffectPositionNormalColor->SetVertexColorEnabled(true); // triggers usage of position/normal/color vertices

		// Input layout
		void const* shaderByteCode;
		size_t byteCodeLength;
		g_pEffectPositionNormalColor->GetVertexShaderBytecode(&shaderByteCode, &byteCodeLength);

		pd3dDevice->CreateInputLayout(VertexPositionNormalColor::InputElements,
			VertexPositionNormalColor::InputElementCount,
			shaderByteCode, byteCodeLength,
			&g_pInputLayoutPositionNormalColor);

		// Primitive batch
		g_pPrimitiveBatchPositionNormalColor = new PrimitiveBatch<VertexPositionNormalColor>(pd3dImmediateContext);
	}

	return S_OK;
}

//--------------------------------------------------------------------------------------
// Release D3D11 resources created in OnD3D11CreateDevice 
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D11DestroyDevice(void* pUserContext)
{
	SAFE_RELEASE(g_pEffect);

	TwDeleteBar(g_pTweakBar);
	g_pTweakBar = nullptr;
	TwTerminate();

	g_pSphere.reset();
	g_pTeapot.reset();
	g_pCube.reset();

	SAFE_DELETE(g_pPrimitiveBatchPositionColor);
	SAFE_RELEASE(g_pInputLayoutPositionColor);
	SAFE_DELETE(g_pEffectPositionColor);

	SAFE_DELETE(g_pPrimitiveBatchPositionNormal);
	SAFE_RELEASE(g_pInputLayoutPositionNormal);
	SAFE_DELETE(g_pEffectPositionNormal);

	SAFE_DELETE(g_pPrimitiveBatchPositionNormalColor);
	SAFE_RELEASE(g_pInputLayoutPositionNormalColor);
	SAFE_DELETE(g_pEffectPositionNormalColor);
}

//--------------------------------------------------------------------------------------
// Create any D3D11 resources that depend on the back buffer
//--------------------------------------------------------------------------------------
HRESULT CALLBACK OnD3D11ResizedSwapChain(ID3D11Device* pd3dDevice, IDXGISwapChain* pSwapChain,
	const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc, void* pUserContext)
{
	// Update camera parameters
	int width = pBackBufferSurfaceDesc->Width;
	int height = pBackBufferSurfaceDesc->Height;
	g_camera.SetWindow(width, height);
	g_camera.SetProjParams(XM_PI / 4.0f, float(width) / float(height), 0.1f, 100.0f);

	// Inform AntTweakBar about back buffer resolution change
	TwWindowSize(pBackBufferSurfaceDesc->Width, pBackBufferSurfaceDesc->Height);

	return S_OK;
}

//--------------------------------------------------------------------------------------
// Release D3D11 resources created in OnD3D11ResizedSwapChain 
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D11ReleasingSwapChain(void* pUserContext)
{
}

//--------------------------------------------------------------------------------------
// Handle key presses
//--------------------------------------------------------------------------------------
void CALLBACK OnKeyboard(UINT nChar, bool bKeyDown, bool bAltDown, void* pUserContext)
{
	HRESULT hr;

	if (bKeyDown)
	{
		switch (nChar)
		{
			// RETURN: toggle fullscreen
		case VK_RETURN:
		{
			if (bAltDown) DXUTToggleFullScreen();
			break;
		}
			// F8: Take screenshot
		case VK_F8:
		{
			// Save current render target as png
			static int nr = 0;
			std::wstringstream ss;
			ss << L"Screenshot" << std::setfill(L'0') << std::setw(4) << nr++ << L".png";

			ID3D11Resource* pTex2D = nullptr;
			DXUTGetD3D11RenderTargetView()->GetResource(&pTex2D);
			SaveWICTextureToFile(DXUTGetD3D11DeviceContext(), pTex2D, GUID_ContainerFormatPng, ss.str().c_str());
			SAFE_RELEASE(pTex2D);

			std::wcout << L"Screenshot written to " << ss.str() << std::endl;
			break;
		}
			// F10: Toggle video recording
		case VK_F10:
		{
			if (!g_pFFmpegVideoRecorder) {
				g_pFFmpegVideoRecorder = new FFmpeg(25, 21, FFmpeg::MODE_INTERPOLATE);
				V(g_pFFmpegVideoRecorder->StartRecording(DXUTGetD3D11Device(), DXUTGetD3D11RenderTargetView(), "output.avi"));
			}
			else {
				g_pFFmpegVideoRecorder->StopRecording();
				SAFE_DELETE(g_pFFmpegVideoRecorder);
			}
		}
		}
	}
}


//--------------------------------------------------------------------------------------
// Handle mouse button presses
//--------------------------------------------------------------------------------------
void CALLBACK OnMouse(bool bLeftButtonDown, bool bRightButtonDown, bool bMiddleButtonDown,
	bool bSideButton1Down, bool bSideButton2Down, int nMouseWheelDelta,
	int xPos, int yPos, void* pUserContext)
{
	// Track mouse movement if left mouse key is pressed
	{
		static int xPosSave = 0, yPosSave = 0;

		if (bLeftButtonDown)
		{
			// Accumulate deltas in g_viMouseDelta
			g_viMouseDelta.x += xPos - xPosSave;
			g_viMouseDelta.y += yPos - yPosSave;
		}

		xPosSave = xPos;
		yPosSave = yPos;
	}
}


//--------------------------------------------------------------------------------------
// Handle messages to the application
//--------------------------------------------------------------------------------------
LRESULT CALLBACK MsgProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam,
	bool* pbNoFurtherProcessing, void* pUserContext)
{
	// Send message to AntTweakbar first
	if (TwEventWin(hWnd, uMsg, wParam, lParam))
	{
		*pbNoFurtherProcessing = true;
		return 0;
	}

	// If message not processed yet, send to camera
	if (g_camera.HandleMessages(hWnd, uMsg, wParam, lParam))
	{
		*pbNoFurtherProcessing = true;
		return 0;
	}

	return 0;
}


//--------------------------------------------------------------------------------------
// Handle updates to the scene
//--------------------------------------------------------------------------------------
void CALLBACK OnFrameMove(double dTime, float fElapsedTime, void* pUserContext)
{
	UpdateWindowTitle(L"Demo");

	// Move camera
	g_camera.FrameMove(fElapsedTime);

	// Update effects with new view + proj transformations
	g_pEffectPositionColor->SetView(g_camera.GetViewMatrix());
	g_pEffectPositionColor->SetProjection(g_camera.GetProjMatrix());

	g_pEffectPositionNormal->SetView(g_camera.GetViewMatrix());
	g_pEffectPositionNormal->SetProjection(g_camera.GetProjMatrix());

	g_pEffectPositionNormalColor->SetView(g_camera.GetViewMatrix());
	g_pEffectPositionNormalColor->SetProjection(g_camera.GetProjMatrix());

	// Apply accumulated mouse deltas to g_vfMovableObjectPos (move along cameras view plane)
	if (g_viMouseDelta.x != 0 || g_viMouseDelta.y != 0)
	{
		// Calcuate camera directions in world space
		XMMATRIX viewInv = XMMatrixInverse(nullptr, g_camera.GetViewMatrix());
		XMVECTOR camRightWorld = XMVector3TransformNormal(g_XMIdentityR0, viewInv);
		XMVECTOR camUpWorld = XMVector3TransformNormal(g_XMIdentityR1, viewInv);

		// Add accumulated mouse deltas to movable object pos
		XMVECTOR vMovableObjectPos = XMLoadFloat3(&g_vfMovableObjectPos);

		float speedScale = 0.001f;
		vMovableObjectPos = XMVectorAdd(vMovableObjectPos, speedScale * (float)g_viMouseDelta.x * camRightWorld);
		vMovableObjectPos = XMVectorAdd(vMovableObjectPos, -speedScale * (float)g_viMouseDelta.y * camUpWorld);

		XMStoreFloat3(&g_vfMovableObjectPos, vMovableObjectPos);

		// Reset accumulated mouse deltas
		g_viMouseDelta = XMINT2(0, 0);
	}

	f_timeAcc += fElapsedTime;
	while (f_timeAcc > g_fTimeStepSize)
	{
		if (g_bMassSpringSystem){
			if (g_iNumSpheres != i_oldNum || g_fSphereSize != f_oldSize)
			{
				InitPoints();
				InitSprings();
				i_oldNum = g_iNumSpheres;
				f_oldSize = g_fSphereSize;
			}
			SetStiffness(g_fStiffness);
			SetMass(g_fMass);
			ApplyPhysikMSS();
		}
		else if (g_bRigidbody)
		{
			//ApplyGravityToRigidbodys();
			ApplyPhysikRBS();
			CollisionDetectionRigidbody();
		}
		else if (g_bBiab)
		{
			if (g_iNumSpheres != i_oldNum || g_fSphereSize != f_oldSize)
			{
				InitBalls();
				i_oldNum = g_iNumSpheres;
				f_oldSize = g_fSphereSize;
			}
			ApplyGravityToRigidbodysToBalls();
			NaiveCollisionDetection();
			UpdateBallPosition();
			//RestrainingPosition();
			//UpdateBallPosition();
		}
		f_timeAcc -= g_fTimeStepSize;
	}
}

//--------------------------------------------------------------------------------------
// Render the scene using the D3D11 device
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D11FrameRender(ID3D11Device* pd3dDevice, ID3D11DeviceContext* pd3dImmediateContext,
	double fTime, float fElapsedTime, void* pUserContext)
{
	HRESULT hr;

	// Clear render target and depth stencil
	float ClearColor[4] = { 0.0f, 0.0f, 0.0f, 0.0f };

	ID3D11RenderTargetView* pRTV = DXUTGetD3D11RenderTargetView();
	ID3D11DepthStencilView* pDSV = DXUTGetD3D11DepthStencilView();
	pd3dImmediateContext->ClearRenderTargetView(pRTV, ClearColor);
	pd3dImmediateContext->ClearDepthStencilView(pDSV, D3D11_CLEAR_DEPTH, 1.0f, 0);

	// Draw floor
	DrawFloor(pd3dImmediateContext);

	// Draw axis box
	DrawBoundingBox(pd3dImmediateContext);

	// Draw speheres
	if (g_bMassSpringSystem) DrawMassSpringSystem(pd3dImmediateContext);

	// Excercise 2: Rigid bodies
	if (g_bRigidbody) DrawRigidBody(pd3dImmediateContext);

	// Excercise 3: Balls in a box
	if (g_bBiab) DrawBiab(pd3dImmediateContext);

	// Draw GUI
	TwDraw();

	if (g_pFFmpegVideoRecorder)
	{
		V(g_pFFmpegVideoRecorder->AddFrame(pd3dImmediateContext, DXUTGetD3D11RenderTargetView()));
	}
}

//--------------------------------------------------------------------------------------
// Initialize everything and go into a render loop
//--------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
#if defined(DEBUG) | defined(_DEBUG)
	// Enable run-time memory check for debug builds.
	// (on program exit, memory leaks are printed to Visual Studio's Output console)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

#ifdef _DEBUG
	std::wcout << L"---- DEBUG BUILD ----\n\n";
#endif

	// Set general DXUT callbacks
	DXUTSetCallbackMsgProc(MsgProc);
	DXUTSetCallbackMouse(OnMouse, true);
	DXUTSetCallbackKeyboard(OnKeyboard);

	DXUTSetCallbackFrameMove(OnFrameMove);
	DXUTSetCallbackDeviceChanging(ModifyDeviceSettings);

	// Set the D3D11 DXUT callbacks
	DXUTSetCallbackD3D11DeviceAcceptable(IsD3D11DeviceAcceptable);
	DXUTSetCallbackD3D11DeviceCreated(OnD3D11CreateDevice);
	DXUTSetCallbackD3D11SwapChainResized(OnD3D11ResizedSwapChain);
	DXUTSetCallbackD3D11FrameRender(OnD3D11FrameRender);
	DXUTSetCallbackD3D11SwapChainReleasing(OnD3D11ReleasingSwapChain);
	DXUTSetCallbackD3D11DeviceDestroyed(OnD3D11DestroyDevice);

	// Init camera
	XMFLOAT3 eye(0.0f, 0.0f, -2.0f);
	XMFLOAT3 lookAt(0.0f, 0.0f, 0.0f);
	g_camera.SetViewParams(XMLoadFloat3(&eye), XMLoadFloat3(&lookAt));
	g_camera.SetButtonMasks(MOUSE_MIDDLE_BUTTON, MOUSE_WHEEL, MOUSE_RIGHT_BUTTON);

	// Init DXUT and create device
	DXUTInit(true, true, NULL); // Parse the command line, show msgboxes on error, no extra command line params
	//DXUTSetIsInGammaCorrectMode( false ); // true by default (SRGB backbuffer), disable to force a RGB backbuffer
	DXUTSetCursorSettings(true, true); // Show the cursor and clip it when in full screen
	DXUTCreateWindow(L"Demo");
	DXUTCreateDevice(D3D_FEATURE_LEVEL_10_0, true, 1280, 960);

	DXUTMainLoop(); // Enter into the DXUT render loop

	DXUTShutdown(); // Shuts down DXUT (includes calls to OnD3D11ReleasingSwapChain() and OnD3D11DestroyDevice())

	return DXUTGetExitCode();
}
