#include "PhysicsObject.h"
#include "PhysicsSystem.h"
#include "../CSC8503Common/Transform.h"
#include "H:\Downloads\8503 Release  2021\Common\Maths.h"

using namespace NCL;
using namespace Maths;
using namespace CSC8503;

PhysicsObject::PhysicsObject(Transform* parentTransform, const CollisionVolume* parentVolume)	{
	transform	= parentTransform;
	volume		= parentVolume;

	inverseMass = 1.0f;
	elasticity	= 0.8f;
	friction	= 0.8f;
}

PhysicsObject::~PhysicsObject()	{

}

void PhysicsObject::ApplyAngularImpulse(const Vector3& force) {
	if (force.Length() > 0) {
		bool a = true;
	}
	angularVelocity += inverseInteriaTensor * force;
}

void PhysicsObject::ApplyLinearImpulse(const Vector3& force) {
	linearVelocity += force * inverseMass;
}

void PhysicsObject::AddForce(const Vector3& addedForce) {
	force += addedForce;
}

void PhysicsObject::AddForceAtPosition(const Vector3& addedForce, const Vector3& position) {
	Vector3 localPos = position - transform->GetPosition();

	force  += addedForce;
	torque += Vector3::Cross(localPos, addedForce);
}

void PhysicsObject::AddTorque(const Vector3& addedTorque) {
	torque += addedTorque;
}

void PhysicsObject::ClearForces() {
	force				= Vector3();
	torque				= Vector3();
}

void PhysicsObject::InitCubeInertia() {
	Vector3 dimensions	= transform->GetScale();

	Vector3 fullWidth = dimensions * 2;

	Vector3 dimsSqr		= fullWidth * fullWidth;

	inverseInertia.x = (12.0f * inverseMass) / (dimsSqr.y + dimsSqr.z);
	inverseInertia.y = (12.0f * inverseMass) / (dimsSqr.x + dimsSqr.z);
	inverseInertia.z = (12.0f * inverseMass) / (dimsSqr.x + dimsSqr.y);
}

void PhysicsObject::InitCapsuleInertia(float radius, float height) {
	inverseInertia.x = 1 / (height * radius * radius * PI *(((height * height)/12) + ((radius * radius) / 4)) + (4 * radius * radius * radius * PI/3) * (((2 * radius * radius) / 5) + ((height * height) / 2) + ((3 * height * radius) / 8)));
	inverseInertia.y = 1 / (height * radius * radius * PI * ((radius * radius) / 2) + (4 * radius * radius * radius * PI / 3) * ((2 * radius * radius) / 5));
	inverseInertia.z = 1 / (height * radius * radius * PI * (((height * height) / 12) + ((radius * radius) / 4)) + (4 * radius * radius * radius * PI / 3) * (((2 * radius * radius) / 5) + ((height * height) / 2) + ((3 * height * radius) / 8)));
}


void PhysicsObject::InitSphereInertia() {
	float radius	= transform->GetScale().GetMaxElement();
	float i			= 2.5f * inverseMass / (radius*radius);

	inverseInertia	= Vector3(i, i, i);
}

void PhysicsObject::InitHollowSphereInertia() {
	float radius = transform->GetScale().GetMaxElement();
	float i = 1.5f * inverseMass / (radius * radius);

	inverseInertia = Vector3(i, i, i);
}

void PhysicsObject::UpdateInertiaTensor() {
	Quaternion q = transform->GetOrientation();
	
	Matrix3 invOrientation	= Matrix3(q.Conjugate());
	Matrix3 orientation		= Matrix3(q);

	inverseInteriaTensor = orientation * Matrix3::Scale(inverseInertia) *invOrientation;
}