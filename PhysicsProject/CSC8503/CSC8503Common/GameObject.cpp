#include "GameObject.h"
#include "CollisionDetection.h"

using namespace NCL::CSC8503;

GameObject::GameObject(string objectName)	{
	name			= objectName;
	worldID			= -1;
	isActive		= true;
	isTrigger		= false;
	deteriorate		= false;
	boundingVolume	= nullptr;
	physicsObject	= nullptr;
	renderObject	= nullptr;
	type = ObjectType::Nothing;
}

GameObject::~GameObject()	{
	delete boundingVolume;
	delete physicsObject;
	delete renderObject;
}

bool GameObject::GetBroadphaseAABB(Vector3&outSize) const {
	if (!boundingVolume) {
		return false;
	}
	outSize = broadphaseAABB;
	return true;
}

void GameObject::UpdateBroadphaseAABB() {
	if (!boundingVolume) {
		return;
	}
	if (boundingVolume->type == VolumeType::AABB) {
		broadphaseAABB = ((AABBVolume&)*boundingVolume).GetHalfDimensions();
	}
	else if (boundingVolume->type == VolumeType::Sphere) {
		float r = ((SphereVolume&)*boundingVolume).GetRadius();
		broadphaseAABB = Vector3(r, r, r);
	}
	else if (boundingVolume->type == VolumeType::OBB) {
		Matrix3 mat = Matrix3(transform.GetOrientation());
		mat = mat.Absolute();
		Vector3 halfSizes = ((OBBVolume&)*boundingVolume).GetHalfDimensions();
		broadphaseAABB = mat * halfSizes;
	}
	else if (boundingVolume->type == VolumeType::Capsule) {
		float halfheight = ((CapsuleVolume&)* boundingVolume).GetHalfHeight();
		float r = ((CapsuleVolume&)* boundingVolume).GetRadius();
		broadphaseAABB = Vector3(halfheight + r, halfheight + r, halfheight + r);
	}
}