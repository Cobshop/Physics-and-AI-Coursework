#include "CollisionDetection.h"
#include "CollisionVolume.h"
#include "AABBVolume.h"
#include "OBBVolume.h"
#include "SphereVolume.h"
#include "../../Common/Vector2.h"
#include "../../Common/Window.h"
#include "../../Common/Maths.h"
#include "Debug.h"

#include <list>

using namespace NCL;

bool CollisionDetection::RayPlaneIntersection(const Ray&r, const Plane&p, RayCollision& collisions) {
	float ln = Vector3::Dot(p.GetNormal(), r.GetDirection());

	if (ln == 0.0f) {
		return false; //direction vectors are perpendicular!
	}
	
	Vector3 planePoint = p.GetPointOnPlane();

	Vector3 pointDir = planePoint - r.GetPosition();

	float d = Vector3::Dot(pointDir, p.GetNormal()) / ln;

	collisions.collidedAt = r.GetPosition() + (r.GetDirection() * d);

	return true;
}

bool CollisionDetection::RayIntersection(const Ray& r,GameObject& object, RayCollision& collision) {
	bool hasCollided = false;

	const Transform& worldTransform = object.GetTransform();
	const CollisionVolume* volume	= object.GetBoundingVolume();

	if (!volume) {
		return false;
	}

	switch (volume->type) {
		case VolumeType::AABB:		hasCollided = RayAABBIntersection(r, worldTransform, (const AABBVolume&)*volume	, collision); break;
		case VolumeType::OBB:		hasCollided = RayOBBIntersection(r, worldTransform, (const OBBVolume&)*volume	, collision); break;
		case VolumeType::Sphere:	hasCollided = RaySphereIntersection(r, worldTransform, (const SphereVolume&)*volume	, collision); break;
		case VolumeType::Capsule:	hasCollided = RayCapsuleIntersection(r, worldTransform, (const CapsuleVolume&)*volume, collision); break;
	}

	return hasCollided;
}

bool CollisionDetection::RayBoxIntersection(const Ray&r, const Vector3& boxPos, const Vector3& boxSize, RayCollision& collision) {
	Vector3 boxMin = boxPos - boxSize;
	Vector3 boxMax = boxPos + boxSize;

	Vector3 rayPos = r.GetPosition();
	Vector3 rayDir = r.GetDirection();

	Vector3 tVals(-1, -1, -1);

	for (int i = 0; i < 3; ++i) { //get best 3 intersections
		if (rayDir[i] > 0) {
			tVals[i] = (boxMin[i] - rayPos[i]) / rayDir[i];
		}
		else if (rayDir[i] < 0) {
			tVals[i] = (boxMax[i] - rayPos[i]) / rayDir[i];
		}
	}
	float bestT = tVals.GetMaxElement();
	if (bestT < 0.0f) {
		return false; //no backwards rays!
	}

	Vector3 intersection = rayPos + (rayDir * bestT);
	const float epsilon = 0.0001f; //an amount of leeway in our calcs
	for (int i = 0; i < 3; ++i) {
	if (intersection[i] + epsilon < boxMin[i] ||
		intersection[i] - epsilon > boxMax[i]) {
			return false; //best intersection doesn’t touch the box!
		}
	}
	collision.collidedAt = intersection;
	collision.rayDistance = bestT;

	return true;
}

bool CollisionDetection::RayAABBIntersection(const Ray&r, const Transform& worldTransform, const AABBVolume& volume, RayCollision& collision) {
	Vector3 boxPos = worldTransform.GetPosition();
	Vector3 boxSize = volume.GetHalfDimensions();
	return RayBoxIntersection(r, boxPos, boxSize, collision);
}

bool CollisionDetection::RayOBBIntersection(const Ray&r, const Transform& worldTransform, const OBBVolume& volume, RayCollision& collision) {
	Quaternion orientation = worldTransform.GetOrientation();
	Vector3 position = worldTransform.GetPosition();
	
	Matrix3 transform = Matrix3(orientation);
	Matrix3 invTransform = Matrix3(orientation.Conjugate());
	
	Vector3 localRayPos = r.GetPosition() - position;
	
	Ray tempRay(invTransform * localRayPos, invTransform * r.GetDirection());
	
	bool collided = RayBoxIntersection(tempRay, Vector3(), volume.GetHalfDimensions(), collision);
	
	if (collided) {
		collision.collidedAt = transform * collision.collidedAt + position;
	
	}
	return collided;
}

bool CollisionDetection::RayCapsuleIntersection(const Ray& r, const Transform& worldTransform, const CapsuleVolume& volume, RayCollision& collision) {
	float Radius = volume.GetRadius();
	Vector3 ExtremeA = worldTransform.GetPosition() + worldTransform.GetOrientation() * Vector3(0, 1, 0) * (volume.GetHalfHeight() - Radius);
	Vector3 ExtremeB = worldTransform.GetPosition() - worldTransform.GetOrientation() * Vector3(0, 1, 0) * (volume.GetHalfHeight() - Radius);
	Vector3 direction = r.GetDirection();

	Vector3 length = ExtremeA - ExtremeB;
	Vector3 rayToExtremeB = r.GetPosition() - ExtremeB;
	float baba = Vector3::Dot(length, length);
	float bard = Vector3::Dot(length, direction);
	float baoa = Vector3::Dot(length, rayToExtremeB);
	float rdoa = Vector3::Dot(direction, rayToExtremeB);
	float oaoa = Vector3::Dot(rayToExtremeB, rayToExtremeB);
	float a = baba - bard * bard;
	float b = baba * rdoa - baoa * bard;
	float c = baba * oaoa - baoa * baoa - Radius * Radius * baba;
	float h = b * b - a * c;
	if (h >= 0.0)
	{
		float t = (-b - sqrt(h)) / a;
		float y = baoa + t * bard;
		// body
		if (y > 0.0 && y < baba) {
			collision.rayDistance = t;
			collision.collidedAt = r.GetPosition() + (direction * collision.rayDistance);
			return true;
		}
		// caps
		Vector3 oc = (y <= 0.0) ? rayToExtremeB : r.GetPosition() - ExtremeA;
		b = Vector3::Dot(direction, oc);
		c = Vector3::Dot(oc, oc) - Radius * Radius;
		h = b * b - c;
		
		if (h > 0.0) {
			collision.rayDistance = -b - sqrt(h);
			collision.collidedAt = r.GetPosition() + (direction * collision.rayDistance);
			return true;
		}
	}
	return false;
}

bool CollisionDetection::RaySphereIntersection(const Ray&r, const Transform& worldTransform, const SphereVolume& volume, RayCollision& collision) {
	Vector3 spherePos = worldTransform.GetPosition();
	float sphereRadius = volume.GetRadius();
	
	//Get the direction between the ray origin and the sphere origin
	Vector3 dir = (spherePos - r.GetPosition());
	
	//Then project the sphere ’s origin onto our ray direction vector
	float sphereProj = Vector3::Dot(dir, r.GetDirection());
	
	if (sphereProj < 0.0f) {
		return false; // point is behind the ray!
	}
	
	//Get closest point on ray line to sphere
	Vector3 point = r.GetPosition() + (r.GetDirection() * sphereProj);
	
	float sphereDist = (point - spherePos).Length();
	
	if (sphereDist > sphereRadius) {
		return false;
	}
	
	float offset = sqrt((sphereRadius * sphereRadius) - (sphereDist * sphereDist));
	
	collision.rayDistance = sphereProj - (offset);
	collision.collidedAt = r.GetPosition() + (r.GetDirection() * collision.rayDistance);
	return true;
}

Matrix4 GenerateInverseView(const Camera &c) {
	float pitch = c.GetPitch();
	float yaw	= c.GetYaw();
	Vector3 position = c.GetPosition();

	Matrix4 iview =
		Matrix4::Translation(position) *
		Matrix4::Rotation(-yaw, Vector3(0, -1, 0)) *
		Matrix4::Rotation(-pitch, Vector3(-1, 0, 0));

	return iview;
}

Vector3 CollisionDetection::Unproject(const Vector3& screenPos, const Camera& cam) {
	Vector2 screenSize = Window::GetWindow()->GetScreenSize();

	float aspect	= screenSize.x / screenSize.y;
	float fov		= cam.GetFieldOfVision();
	float nearPlane = cam.GetNearPlane();
	float farPlane  = cam.GetFarPlane();

	//Create our inverted matrix! Note how that to get a correct inverse matrix,
	//the order of matrices used to form it are inverted, too.
	Matrix4 invVP = GenerateInverseView(cam) * GenerateInverseProjection(aspect, fov, nearPlane, farPlane);

	//Our mouse position x and y values are in 0 to screen dimensions range,
	//so we need to turn them into the -1 to 1 axis range of clip space.
	//We can do that by dividing the mouse values by the width and height of the
	//screen (giving us a range of 0.0 to 1.0), multiplying by 2 (0.0 to 2.0)
	//and then subtracting 1 (-1.0 to 1.0).
	Vector4 clipSpace = Vector4(
		(screenPos.x / (float)screenSize.x) * 2.0f - 1.0f,
		(screenPos.y / (float)screenSize.y) * 2.0f - 1.0f,
		(screenPos.z),
		1.0f
	);

	//Then, we multiply our clipspace coordinate by our inverted matrix
	Vector4 transformed = invVP * clipSpace;

	//our transformed w coordinate is now the 'inverse' perspective divide, so
	//we can reconstruct the final world space by dividing x,y,and z by w.
	return Vector3(transformed.x / transformed.w, transformed.y / transformed.w, transformed.z / transformed.w);
}

Ray CollisionDetection::BuildRayFromMouse(const Camera& cam) {
	Vector2 screenMouse = Window::GetMouse()->GetAbsolutePosition();
	Vector2 screenSize	= Window::GetWindow()->GetScreenSize();

	//We remove the y axis mouse position from height as OpenGL is 'upside down',
	//and thinks the bottom left is the origin, instead of the top left!
	Vector3 nearPos = Vector3(screenMouse.x,
		screenSize.y - screenMouse.y,
		-0.99999f
	);

	//We also don't use exactly 1.0 (the normalised 'end' of the far plane) as this
	//causes the unproject function to go a bit weird. 
	Vector3 farPos = Vector3(screenMouse.x,
		screenSize.y - screenMouse.y,
		0.99999f
	);

	Vector3 a = Unproject(nearPos, cam);
	Vector3 b = Unproject(farPos, cam);
	Vector3 c = b - a;

	c.Normalise();

	//std::cout << "Ray Direction:" << c << std::endl;

	return Ray(cam.GetPosition(), c);
}

//http://bookofhook.com/mousepick.pdf
Matrix4 CollisionDetection::GenerateInverseProjection(float aspect, float fov, float nearPlane, float farPlane) {
	Matrix4 m;

	float t = tan(fov*PI_OVER_360);

	float neg_depth = nearPlane - farPlane;

	const float h = 1.0f / t;

	float c = (farPlane + nearPlane) / neg_depth;
	float e = -1.0f;
	float d = 2.0f*(nearPlane*farPlane) / neg_depth;

	m.array[0]  = aspect / h;
	m.array[5]  = tan(fov*PI_OVER_360);

	m.array[10] = 0.0f;
	m.array[11] = 1.0f / d;

	m.array[14] = 1.0f / e;

	m.array[15] = -c / (d*e);

	return m;
}

/*
And here's how we generate an inverse view matrix. It's pretty much
an exact inversion of the BuildViewMatrix function of the Camera class!
*/
Matrix4 CollisionDetection::GenerateInverseView(const Camera &c) {
	float pitch = c.GetPitch();
	float yaw	= c.GetYaw();
	Vector3 position = c.GetPosition();

	Matrix4 iview =
Matrix4::Translation(position) *
Matrix4::Rotation(yaw, Vector3(0, 1, 0)) *
Matrix4::Rotation(pitch, Vector3(1, 0, 0));

return iview;
}


/*
If you've read through the Deferred Rendering tutorial you should have a pretty
good idea what this function does. It takes a 2D position, such as the mouse
position, and 'unprojects' it, to generate a 3D world space position for it.

Just as we turn a world space position into a clip space position by multiplying
it by the model, view, and projection matrices, we can turn a clip space
position back to a 3D position by multiply it by the INVERSE of the
view projection matrix (the model matrix has already been assumed to have
'transformed' the 2D point). As has been mentioned a few times, inverting a
matrix is not a nice operation, either to understand or code. But! We can cheat
the inversion process again, just like we do when we create a view matrix using
the camera.

So, to form the inverted matrix, we need the aspect and fov used to create the
projection matrix of our scene, and the camera used to form the view matrix.

*/
Vector3	CollisionDetection::UnprojectScreenPosition(Vector3 position, float aspect, float fov, const Camera &c) {
	//Create our inverted matrix! Note how that to get a correct inverse matrix,
	//the order of matrices used to form it are inverted, too.
	Matrix4 invVP = GenerateInverseView(c) * GenerateInverseProjection(aspect, fov, c.GetNearPlane(), c.GetFarPlane());

	Vector2 screenSize = Window::GetWindow()->GetScreenSize();

	//Our mouse position x and y values are in 0 to screen dimensions range,
	//so we need to turn them into the -1 to 1 axis range of clip space.
	//We can do that by dividing the mouse values by the width and height of the
	//screen (giving us a range of 0.0 to 1.0), multiplying by 2 (0.0 to 2.0)
	//and then subtracting 1 (-1.0 to 1.0).
	Vector4 clipSpace = Vector4(
		(position.x / (float)screenSize.x) * 2.0f - 1.0f,
		(position.y / (float)screenSize.y) * 2.0f - 1.0f,
		(position.z) - 1.0f,
		1.0f
	);

	//Then, we multiply our clipspace coordinate by our inverted matrix
	Vector4 transformed = invVP * clipSpace;

	//our transformed w coordinate is now the 'inverse' perspective divide, so
	//we can reconstruct the final world space by dividing x,y,and z by w.
	return Vector3(transformed.x / transformed.w, transformed.y / transformed.w, transformed.z / transformed.w);
}

bool CollisionDetection::ObjectIntersection(GameObject* a, GameObject* b, CollisionInfo& collisionInfo) {
	const CollisionVolume* volA = a->GetBoundingVolume();
	const CollisionVolume* volB = b->GetBoundingVolume();

	if (!volA || !volB) {
		return false;
	}

	if (a->GetPhysicsObject()->GetInverseMass() == 0 && b->GetPhysicsObject()->GetInverseMass() == 0) {
		return false;
	}

	collisionInfo.a = a;
	collisionInfo.b = b;

	Transform& transformA = a->GetTransform();
	Transform& transformB = b->GetTransform();

	VolumeType pairType = (VolumeType)((int)volA->type | (int)volB->type);

	if (pairType == VolumeType::AABB) {
		return AABBIntersection((AABBVolume&)*volA, transformA, (AABBVolume&)*volB, transformB, collisionInfo);
	}

	if (pairType == VolumeType::Sphere) {
		return SphereIntersection((SphereVolume&)*volA, transformA, (SphereVolume&)*volB, transformB, collisionInfo);
	}

	if (pairType == VolumeType::OBB) {
		return OBBIntersection((OBBVolume&)*volA, transformA, (OBBVolume&)*volB, transformB, collisionInfo);
	}

	if (pairType == VolumeType::Capsule) {
		return CapsuleCapsuleIntersection((CapsuleVolume&)* volA, transformA, (CapsuleVolume&)* volB, transformB, collisionInfo);
	}

	if (volA->type == VolumeType::AABB && volB->type == VolumeType::Sphere) {
		return AABBSphereIntersection((AABBVolume&)*volA, transformA, (SphereVolume&)*volB, transformB, collisionInfo, false);
	}
	if (volA->type == VolumeType::Sphere && volB->type == VolumeType::AABB) {
		collisionInfo.a = b;
		collisionInfo.b = a;
		return AABBSphereIntersection((AABBVolume&)*volB, transformB, (SphereVolume&)*volA, transformA, collisionInfo, false);
	}

	if (volA->type == VolumeType::OBB && volB->type == VolumeType::Sphere) {
		return OBBSphereIntersection((OBBVolume&)* volA, transformA, (SphereVolume&)* volB, transformB, collisionInfo, false);
	}
	if (volA->type == VolumeType::Sphere && volB->type == VolumeType::OBB) {
		collisionInfo.a = b;
		collisionInfo.b = a;
		return OBBSphereIntersection((OBBVolume&)* volB, transformB, (SphereVolume&)* volA, transformA, collisionInfo, false);
	}

	if (volA->type == VolumeType::Capsule && volB->type == VolumeType::Sphere) {
		return SphereCapsuleIntersection((CapsuleVolume&)*volA, transformA, (SphereVolume&)*volB, transformB, collisionInfo);
	}
	if (volA->type == VolumeType::Sphere && volB->type == VolumeType::Capsule) {
		collisionInfo.a = b;
		collisionInfo.b = a;
		return SphereCapsuleIntersection((CapsuleVolume&)*volB, transformB, (SphereVolume&)*volA, transformA, collisionInfo);
	}

	if (volA->type == VolumeType::Capsule && volB->type == VolumeType::AABB) {
		return AABBCapsuleIntersection((CapsuleVolume&)* volA, transformA, (AABBVolume&)* volB, transformB, collisionInfo);
	}
	if (volA->type == VolumeType::AABB && volB->type == VolumeType::Capsule) {
		collisionInfo.a = b;
		collisionInfo.b = a;
		return AABBCapsuleIntersection((CapsuleVolume&)* volB, transformB, (AABBVolume&)* volA, transformA, collisionInfo);
	}

	if (volA->type == VolumeType::Capsule && volB->type == VolumeType::OBB) {
		return OBBCapsuleIntersection((CapsuleVolume&)* volA, transformA, (OBBVolume&)* volB, transformB, collisionInfo);
	}
	if (volA->type == VolumeType::OBB && volB->type == VolumeType::Capsule) {
		collisionInfo.a = b;
		collisionInfo.b = a;
		return OBBCapsuleIntersection((CapsuleVolume&)* volB, transformB, (OBBVolume&)* volA, transformA, collisionInfo);
	}

	if (volA->type == VolumeType::AABB && volB->type == VolumeType::OBB) {
		return AABBOBBIntersection((AABBVolume&)* volA, transformA, (OBBVolume&)* volB, transformB, collisionInfo);
	}
	if (volA->type == VolumeType::OBB && volB->type == VolumeType::AABB) {
		collisionInfo.a = b;
		collisionInfo.b = a;
		return AABBOBBIntersection((AABBVolume&)* volB, transformB, (OBBVolume&)* volA, transformA, collisionInfo);
	}

	return false;
}

bool CollisionDetection::AABBTest(const Vector3& posA, const Vector3& posB, const Vector3& halfSizeA, const Vector3& halfSizeB) {

	Vector3 delta = posB - posA;
	Vector3 totalSize = halfSizeA + halfSizeB;
	
	if (abs(delta.x) < totalSize.x && abs(delta.y) < totalSize.y && abs(delta.z) < totalSize.z) {
		return true;	
	}

	return false;
}

//AABB/AABB Collisions
bool CollisionDetection::AABBIntersection(const AABBVolume& volumeA, const Transform& worldTransformA,
	const AABBVolume& volumeB, const Transform& worldTransformB, CollisionInfo& collisionInfo) {

	Vector3 boxAPos = worldTransformA.GetPosition();
	Vector3 boxBPos = worldTransformB.GetPosition();
	
	Vector3 boxASize = volumeA.GetHalfDimensions();
	Vector3 boxBSize = volumeB.GetHalfDimensions();
	
	bool overlap = AABBTest(boxAPos, boxBPos, boxASize, boxBSize);

	if (overlap) {
		static const Vector3 faces[6] ={Vector3(-1, 0, 0), Vector3(1, 0, 0), Vector3(0, -1, 0), Vector3(0, 1, 0), Vector3(0, 0, -1), Vector3(0, 0, 1),};
		
		Vector3 maxA = boxAPos + boxASize;
		Vector3 minA = boxAPos - boxASize;
		
		Vector3 maxB = boxBPos + boxBSize;
		Vector3 minB = boxBPos - boxBSize;
		
		float distances[6] = {
			(maxB.x - minA.x),// distance of box ’b’ to ’left’ of ’a’.
			(maxA.x - minB.x),// distance of box ’b’ to ’right’ of ’a’.
			(maxB.y - minA.y),// distance of box ’b’ to ’bottom ’ of ’a’.
			(maxA.y - minB.y),// distance of box ’b’ to ’top’ of ’a’.
			(maxB.z - minA.z),// distance of box ’b’ to ’far’ of ’a’.
			(maxA.z - minB.z) // distance of box ’b’ to ’near’ of ’a’.
		};
		float penetration = FLT_MAX;
		Vector3 bestAxis;
		
		for (int i = 0; i < 6; i++) {
			if (distances[i] < penetration) {
				penetration = distances[i];
				bestAxis = faces[i];
			}
		}
		collisionInfo.AddContactPoint(Vector3(), Vector3(), bestAxis, penetration);
		return true;
	}

	return false;
}

//Sphere / Sphere Collision
bool CollisionDetection::SphereIntersection(const SphereVolume& volumeA, const Transform& worldTransformA,
	const SphereVolume& volumeB, const Transform& worldTransformB, CollisionInfo& collisionInfo) {

	float radii = volumeA.GetRadius() + volumeB.GetRadius();
	Vector3 delta = worldTransformB.GetPosition() - worldTransformA.GetPosition();

	float deltaLength = delta.Length();

	if (deltaLength < radii) {
		float penetration = (radii - deltaLength);
		Vector3 normal = delta.Normalised();
		Vector3 localA = normal * volumeA.GetRadius();
		Vector3 localB = -normal * volumeB.GetRadius();

		collisionInfo.AddContactPoint(localA, localB, normal, penetration);

		return true;//we’re colliding!
	}


	return false;
}

//AABB - Sphere Collision
bool CollisionDetection::AABBSphereIntersection(const AABBVolume& volumeA, const Transform& worldTransformA,
	const SphereVolume& volumeB, const Transform& worldTransformB, CollisionInfo& collisionInfo, bool capsule) {

	Vector3 boxSize = volumeA.GetHalfDimensions();
	
	Vector3 delta = worldTransformB.GetPosition() - worldTransformA.GetPosition();
	
	Vector3 closestPointOnBox = Maths::Clamp(delta, -boxSize, boxSize);
	
	Vector3 localPoint = delta - closestPointOnBox;
	float distance = localPoint.Length();
	
	if (distance < volumeB.GetRadius()) {//yes , we’re colliding!
		Vector3 collisionNormal = localPoint.Normalised();
		float penetration = (volumeB.GetRadius() - distance);
		
		Vector3 localA = Vector3();
		Vector3 localB = -collisionNormal * volumeB.GetRadius();
		
		if(capsule)
			collisionInfo.AddContactPoint(localA, localB, -collisionNormal, penetration);
		else
			collisionInfo.AddContactPoint(localA, localB, collisionNormal, penetration);

		return true;
		
	}

	return false;
}

//OBB - Sphere Collision
bool CollisionDetection::OBBSphereIntersection(const OBBVolume& volumeA, const Transform& worldTransformA,
	const SphereVolume& volumeB, const Transform& worldTransformB, CollisionInfo& collisionInfo, bool capsule) {

	Quaternion orientation = worldTransformA.GetOrientation();
	Vector3 position = worldTransformA.GetPosition();

	Matrix3 transform = Matrix3(orientation);
	Matrix3 invTransform = Matrix3(orientation.Conjugate());

	Vector3 localSpherePos = worldTransformB.GetPosition() - position;

	Vector3 tempPosition = invTransform * localSpherePos;

	Vector3 boxSize = volumeA.GetHalfDimensions();

	Vector3 closestPointOnBox = Maths::Clamp(tempPosition, -boxSize, boxSize);

	Vector3 localPoint = tempPosition - closestPointOnBox;
	float distance = localPoint.Length();

	if (distance < volumeB.GetRadius()) {//yes , we’re colliding!
		Vector3 tClosestPointOnBox = transform * closestPointOnBox;
		Vector3 localPointer = localSpherePos - tClosestPointOnBox;
		Vector3 collisionNormal = localPointer.Normalised();
		float penetration = (volumeB.GetRadius() - distance);

		Vector3 localA = tClosestPointOnBox;
		Vector3 localB = -collisionNormal * volumeB.GetRadius();

		if(capsule)
			collisionInfo.AddContactPoint(localA, localB, -collisionNormal, penetration);
		else
			collisionInfo.AddContactPoint(localA, localB, collisionNormal, penetration);

		return true;
	}

	return false;
}

bool CollisionDetection::OBBIntersection(
	const OBBVolume& volumeA, const Transform& worldTransformA,
	const OBBVolume& volumeB, const Transform& worldTransformB, CollisionInfo& collisionInfo) {

	Quaternion orientationA = worldTransformA.GetOrientation();
	Matrix3 invTransformA = Matrix3(orientationA.Conjugate());
	Vector3 positionA = worldTransformA.GetPosition();
	Matrix3 transformA = Matrix3(orientationA);

	Quaternion orientationB = worldTransformB.GetOrientation();
	Matrix3 invTransformB = Matrix3(orientationB.Conjugate());
	Vector3 positionB = worldTransformB.GetPosition();
	Matrix3 transformB = Matrix3(orientationB);

	Vector3 boxAPoints[8];
	boxAPoints[0] = positionA + transformA * Vector3(-volumeA.GetHalfDimensions().x, -volumeA.GetHalfDimensions().y, -volumeA.GetHalfDimensions().z);
	boxAPoints[1] = positionA + transformA * Vector3(volumeA.GetHalfDimensions().x, -volumeA.GetHalfDimensions().y, -volumeA.GetHalfDimensions().z);
	boxAPoints[2] = positionA + transformA * Vector3(-volumeA.GetHalfDimensions().x, volumeA.GetHalfDimensions().y, -volumeA.GetHalfDimensions().z);
	boxAPoints[3] = positionA + transformA * Vector3(-volumeA.GetHalfDimensions().x, -volumeA.GetHalfDimensions().y, volumeA.GetHalfDimensions().z);
	boxAPoints[4] = positionA + transformA * Vector3(volumeA.GetHalfDimensions().x, -volumeA.GetHalfDimensions().y, volumeA.GetHalfDimensions().z);
	boxAPoints[5] = positionA + transformA * Vector3(volumeA.GetHalfDimensions().x, volumeA.GetHalfDimensions().y, -volumeA.GetHalfDimensions().z);
	boxAPoints[6] = positionA + transformA * Vector3(-volumeA.GetHalfDimensions().x, volumeA.GetHalfDimensions().y, volumeA.GetHalfDimensions().z);
	boxAPoints[7] = positionA + transformA * Vector3(volumeA.GetHalfDimensions().x, volumeA.GetHalfDimensions().y, volumeA.GetHalfDimensions().z);


	Vector3 boxBPoints[8];
	boxBPoints[0] = positionB + transformB * Vector3(-volumeB.GetHalfDimensions().x, -volumeB.GetHalfDimensions().y, -volumeB.GetHalfDimensions().z);
	boxBPoints[1] = positionB + transformB * Vector3(volumeB.GetHalfDimensions().x, -volumeB.GetHalfDimensions().y, -volumeB.GetHalfDimensions().z);
	boxBPoints[2] = positionB + transformB * Vector3(-volumeB.GetHalfDimensions().x, volumeB.GetHalfDimensions().y, -volumeB.GetHalfDimensions().z);
	boxBPoints[3] = positionB + transformB * Vector3(-volumeB.GetHalfDimensions().x, -volumeB.GetHalfDimensions().y, volumeB.GetHalfDimensions().z);
	boxBPoints[4] = positionB + transformB * Vector3(volumeB.GetHalfDimensions().x, -volumeB.GetHalfDimensions().y, volumeB.GetHalfDimensions().z);
	boxBPoints[5] = positionB + transformB * Vector3(volumeB.GetHalfDimensions().x, volumeB.GetHalfDimensions().y, -volumeB.GetHalfDimensions().z);
	boxBPoints[6] = positionB + transformB * Vector3(-volumeB.GetHalfDimensions().x, volumeB.GetHalfDimensions().y, volumeB.GetHalfDimensions().z);
	boxBPoints[7] = positionB + transformB * Vector3(-volumeB.GetHalfDimensions().x, volumeB.GetHalfDimensions().y, volumeB.GetHalfDimensions().z);

	Vector3 boxAxes[15];
	boxAxes[0] = transformB * Vector3(1, 0, 0);
	boxAxes[1] = transformB * Vector3(0, 1, 0);
	boxAxes[2] = transformB * Vector3(0, 0, 1);
	boxAxes[3] = transformA * Vector3(1, 0, 0);
	boxAxes[4] = transformA * Vector3(0, 1, 0);
	boxAxes[5] = transformA * Vector3(0, 0, 1);
	boxAxes[6] = Vector3::Cross(boxAxes[0], boxAxes[3]);
	boxAxes[7] = Vector3::Cross(boxAxes[0], boxAxes[4]);
	boxAxes[8] = Vector3::Cross(boxAxes[0], boxAxes[5]);
	boxAxes[9] = Vector3::Cross(boxAxes[1], boxAxes[3]);
	boxAxes[10] = Vector3::Cross(boxAxes[1], boxAxes[4]);
	boxAxes[11] = Vector3::Cross(boxAxes[1], boxAxes[5]);
	boxAxes[12] = Vector3::Cross(boxAxes[2], boxAxes[3]);
	boxAxes[13] = Vector3::Cross(boxAxes[2], boxAxes[4]);
	boxAxes[14] = Vector3::Cross(boxAxes[2], boxAxes[5]);

	Axis Axes[30];
	float penetration = FLT_MAX;
	vector<Vector3> localA;
	vector<Vector3> localB;
	vector<Vector3> normals;
	vector<float> penetrations;
	Vector3 normal;
	bool overlap = true;

	for (int i = 0; i < 15; i++) {
		Axes[i] = ProjectOntoAxis(boxAPoints, boxAxes[i].Normalised());
		Axes[i + 15] = ProjectOntoAxis(boxBPoints, boxAxes[i].Normalised());
	}
	for (int i = 0; i < 15; i++) {
		if (Axes[i].maxProjection < Axes[i + 15].minProjection || Axes[i + 15].maxProjection < Axes[i].minProjection)
			overlap = false;
		else {
			if (Axes[i + 15].maxProjection - Axes[i].minProjection < penetration && Axes[i + 15].maxProjection - Axes[i].minProjection != 0) {
				penetration = Axes[i + 15].maxProjection - Axes[i].minProjection;
				normal = boxAxes[i].Normalised();
				localA.push_back(Axes[i].minPoint);
				localB.push_back(Axes[i + 15].maxPoint);
			}
			if (Axes[i].maxProjection - Axes[i + 15].minProjection < penetration && Axes[i].maxProjection - Axes[i + 15].minProjection != 0) {
				penetration = Axes[i].maxProjection - Axes[i + 15].minProjection;
				normal = boxAxes[i].Normalised();
				localA.push_back(Axes[i].maxPoint);
				localB.push_back(Axes[i + 15].minPoint);
			}
		}
	}

	if (!overlap)
	{
		return false;
	}
	else
	{
		Vector3 tempPosition = invTransformA * (positionB - positionA);
		Vector3 closestPointOnBoxA = Maths::Clamp(tempPosition, -volumeA.GetHalfDimensions(), volumeA.GetHalfDimensions());
		Vector3 localPointA = tempPosition - closestPointOnBoxA;
		float distanceA = localPointA.Length();
		Vector3 tClosestPointOnBoxA = transformA * closestPointOnBoxA;

		tempPosition = invTransformB * (positionA - positionB);
		Vector3 closestPointOnBoxB = Maths::Clamp(tempPosition, -volumeB.GetHalfDimensions(), volumeB.GetHalfDimensions());
		Vector3 tClosestPointOnBoxB = transformB * closestPointOnBoxB;
		Vector3 localPointB = tempPosition - closestPointOnBoxB;
		float distanceB = localPointB.Length();
		//normal = (tClosestPointOnBoxB - tClosestPointOnBoxA).Normalised();

		if (localA.size() < 3)
			collisionInfo.AddContactPoint(tClosestPointOnBoxA, tClosestPointOnBoxB, normal, penetration);
		else
			collisionInfo.AddContactPoint(Vector3(), Vector3(), normal, penetration);

		return true;
	}
}

bool CollisionDetection::AABBOBBIntersection(
	const AABBVolume& volumeA, const Transform& worldTransformA,
	const OBBVolume& volumeB, const Transform& worldTransformB, CollisionInfo& collisionInfo) {

	Quaternion orientationA = worldTransformA.GetOrientation();
	Matrix3 invTransformA = Matrix3(orientationA.Conjugate());
	Vector3 positionA = worldTransformA.GetPosition();
	Matrix3 transformA = Matrix3(orientationA);

	Quaternion orientationB = worldTransformB.GetOrientation();
	Matrix3 invTransformB = Matrix3(orientationB.Conjugate());
	Vector3 positionB = worldTransformB.GetPosition();
	Matrix3 transformB = Matrix3(orientationB);

	Vector3 boxAPoints[8];
	boxAPoints[0] = positionA + transformA * Vector3(-volumeA.GetHalfDimensions().x, -volumeA.GetHalfDimensions().y, -volumeA.GetHalfDimensions().z);
	boxAPoints[1] = positionA + transformA * Vector3(volumeA.GetHalfDimensions().x, -volumeA.GetHalfDimensions().y, -volumeA.GetHalfDimensions().z);
	boxAPoints[2] = positionA + transformA * Vector3(-volumeA.GetHalfDimensions().x, volumeA.GetHalfDimensions().y, -volumeA.GetHalfDimensions().z);
	boxAPoints[3] = positionA + transformA * Vector3(-volumeA.GetHalfDimensions().x, -volumeA.GetHalfDimensions().y, volumeA.GetHalfDimensions().z);
	boxAPoints[4] = positionA + transformA * Vector3(volumeA.GetHalfDimensions().x, -volumeA.GetHalfDimensions().y, volumeA.GetHalfDimensions().z);
	boxAPoints[5] = positionA + transformA * Vector3(volumeA.GetHalfDimensions().x, volumeA.GetHalfDimensions().y, -volumeA.GetHalfDimensions().z);
	boxAPoints[6] = positionA + transformA * Vector3(-volumeA.GetHalfDimensions().x, volumeA.GetHalfDimensions().y, volumeA.GetHalfDimensions().z);
	boxAPoints[7] = positionA + transformA * Vector3(volumeA.GetHalfDimensions().x, volumeA.GetHalfDimensions().y, volumeA.GetHalfDimensions().z);


	Vector3 boxBPoints[8];
	boxBPoints[0] = positionB + transformB * Vector3(-volumeB.GetHalfDimensions().x, -volumeB.GetHalfDimensions().y, -volumeB.GetHalfDimensions().z);
	boxBPoints[1] = positionB + transformB * Vector3(volumeB.GetHalfDimensions().x, -volumeB.GetHalfDimensions().y, -volumeB.GetHalfDimensions().z);
	boxBPoints[2] = positionB + transformB * Vector3(-volumeB.GetHalfDimensions().x, volumeB.GetHalfDimensions().y, -volumeB.GetHalfDimensions().z);
	boxBPoints[3] = positionB + transformB * Vector3(-volumeB.GetHalfDimensions().x, -volumeB.GetHalfDimensions().y, volumeB.GetHalfDimensions().z);
	boxBPoints[4] = positionB + transformB * Vector3(volumeB.GetHalfDimensions().x, -volumeB.GetHalfDimensions().y, volumeB.GetHalfDimensions().z);
	boxBPoints[5] = positionB + transformB * Vector3(volumeB.GetHalfDimensions().x, volumeB.GetHalfDimensions().y, -volumeB.GetHalfDimensions().z);
	boxBPoints[6] = positionB + transformB * Vector3(-volumeB.GetHalfDimensions().x, volumeB.GetHalfDimensions().y, volumeB.GetHalfDimensions().z);
	boxBPoints[7] = positionB + transformB * Vector3(-volumeB.GetHalfDimensions().x, volumeB.GetHalfDimensions().y, volumeB.GetHalfDimensions().z);

	Vector3 boxAxes[15];
	boxAxes[0] = transformB * Vector3(1, 0, 0);
	boxAxes[1] = transformB * Vector3(0, 1, 0);
	boxAxes[2] = transformB * Vector3(0, 0, 1);
	boxAxes[3] = transformA * Vector3(1, 0, 0);
	boxAxes[4] = transformA * Vector3(0, 1, 0);
	boxAxes[5] = transformA * Vector3(0, 0, 1);
	boxAxes[6] = Vector3::Cross(boxAxes[0], boxAxes[3]);
	boxAxes[7] = Vector3::Cross(boxAxes[0], boxAxes[4]);
	boxAxes[8] = Vector3::Cross(boxAxes[0], boxAxes[5]);
	boxAxes[9] = Vector3::Cross(boxAxes[1], boxAxes[3]);
	boxAxes[10] = Vector3::Cross(boxAxes[1], boxAxes[4]);
	boxAxes[11] = Vector3::Cross(boxAxes[1], boxAxes[5]);
	boxAxes[12] = Vector3::Cross(boxAxes[2], boxAxes[3]);
	boxAxes[13] = Vector3::Cross(boxAxes[2], boxAxes[4]);
	boxAxes[14] = Vector3::Cross(boxAxes[2], boxAxes[5]);

	Axis Axes[30];
	float penetration = FLT_MAX;
	vector<Vector3> localA;
	vector<Vector3> localB;
	vector<Vector3> normals;
	vector<float> penetrations;
	Vector3 normal;
	bool overlap = true;

	for (int i = 0; i < 15; i++) {
		Axes[i] = ProjectOntoAxis(boxAPoints, boxAxes[i].Normalised());
		Axes[i + 15] = ProjectOntoAxis(boxBPoints, boxAxes[i].Normalised());
	}
	for (int i = 0; i < 15; i++) {
		if (Axes[i].maxProjection < Axes[i + 15].minProjection || Axes[i + 15].maxProjection < Axes[i].minProjection)
			overlap = false;
		else {
			if (Axes[i + 15].maxProjection - Axes[i].minProjection < penetration && Axes[i + 15].maxProjection - Axes[i].minProjection != 0) {
				penetration = Axes[i + 15].maxProjection - Axes[i].minProjection;
				normal = boxAxes[i].Normalised();
				localA.push_back(Axes[i].minPoint);
				localB.push_back(Axes[i + 15].maxPoint);
			}
			if (Axes[i].maxProjection - Axes[i + 15].minProjection < penetration && Axes[i].maxProjection - Axes[i + 15].minProjection != 0) {
				penetration = Axes[i].maxProjection - Axes[i + 15].minProjection;
				normal = boxAxes[i].Normalised();
				localA.push_back(Axes[i].maxPoint);
				localB.push_back(Axes[i + 15].minPoint);
			}
		}
	}

	if (!overlap)
	{
		return false;
	}
	else
	{
		Vector3 tempPosition = invTransformA * (positionB - positionA);
		Vector3 closestPointOnBoxA = Maths::Clamp(tempPosition, -volumeA.GetHalfDimensions(), volumeA.GetHalfDimensions());
		Vector3 localPointA = tempPosition - closestPointOnBoxA;
		float distanceA = localPointA.Length();
		Vector3 tClosestPointOnBoxA = transformA * closestPointOnBoxA;

		tempPosition = invTransformB * (positionA - positionB);
		Vector3 closestPointOnBoxB = Maths::Clamp(tempPosition, -volumeB.GetHalfDimensions(), volumeB.GetHalfDimensions());
		Vector3 tClosestPointOnBoxB = transformB * closestPointOnBoxB;
		Vector3 localPointB = tempPosition - closestPointOnBoxB;
		float distanceB = localPointB.Length();
		//normal = (tClosestPointOnBoxB - tClosestPointOnBoxA).Normalised();

		if(localA.size() < 3)
			collisionInfo.AddContactPoint(tClosestPointOnBoxA, tClosestPointOnBoxB, normal, penetration);
		else
			collisionInfo.AddContactPoint(Vector3(), Vector3(), normal, penetration);

		return true;
	}
}

CollisionDetection::Axis CollisionDetection::ProjectOntoAxis(const Vector3* boxPoints, const Vector3 normal)
{
	float minProjection = Vector3::Dot(boxPoints[0], normal);
	float maxProjection = minProjection;
	Axis temp;

	for (int i = 0; i < 8; i++)
	{
		float currentProjection = Vector3::Dot(boxPoints[i], normal);

		if (minProjection > currentProjection) {
			minProjection = currentProjection;
			temp.minPoint = boxPoints[i];
		}
		if (currentProjection > maxProjection) {
			maxProjection = currentProjection;
			temp.maxPoint = boxPoints[i];
		}
	}
	temp.maxProjection = maxProjection;
	temp.minProjection = minProjection;
	return temp;
}

bool CollisionDetection::SphereCapsuleIntersection(
	const CapsuleVolume& volumeA, const Transform& worldTransformA,
	const SphereVolume& volumeB, const Transform& worldTransformB, CollisionInfo& collisionInfo) {

	float RadiusA = volumeA.GetRadius();
	Vector3 ExtremeA = worldTransformA.GetPosition() + worldTransformA.GetOrientation() * Vector3(0, 1, 0) * (volumeA.GetHalfHeight() - RadiusA);
	Vector3 ExtremeB = worldTransformA.GetPosition() - worldTransformA.GetOrientation() * Vector3(0, 1, 0) * (volumeA.GetHalfHeight() - RadiusA);
	Vector3 CentreVector = ExtremeB - ExtremeA;

	float p = Vector3::Dot(worldTransformB.GetPosition() - ExtremeA, CentreVector) / Vector3::Dot(CentreVector, CentreVector);
	float distance = Maths::Clamp(p, 0.0f, 1.0f);
	Vector3 ClosestPoint = ExtremeA + (CentreVector * distance);

	SphereVolume ClosestSphere = SphereVolume(RadiusA);
	Transform Transform = worldTransformA;
	Transform.SetPosition(ClosestPoint);

	return SphereIntersection(ClosestSphere, Transform, volumeB, worldTransformB, collisionInfo);
}

bool CollisionDetection::AABBCapsuleIntersection(
	const CapsuleVolume& volumeA, const Transform& worldTransformA,
	const AABBVolume& volumeB, const Transform& worldTransformB, CollisionInfo& collisionInfo) {

	Vector3 position = worldTransformB.GetPosition();
	Vector3 boxSize = volumeB.GetHalfDimensions();
	float RadiusA = volumeA.GetRadius();
	Vector3 ExtremeA = worldTransformA.GetPosition() + worldTransformA.GetOrientation() * Vector3(0, 1, 0) * (volumeA.GetHalfHeight() - RadiusA);
	Vector3 ExtremeB = worldTransformA.GetPosition() - worldTransformA.GetOrientation() * Vector3(0, 1, 0) * (volumeA.GetHalfHeight() - RadiusA);
	Vector3 CentreVector = ExtremeB - ExtremeA;

	Ray CapsuleRay(worldTransformA.GetPosition() - CentreVector * 5000, CentreVector.Normalised());
	vector<Plane> Planes;
	vector<Vector3> normals = {Vector3(1,0,0), Vector3(0,1,0), Vector3(0,0,1)};
	bool collided = false;

	for (int i = 0; i < normals.size(); i++) {
		Planes.push_back(Plane(normals[i], Vector3::Dot(-normals[i], worldTransformB.GetPosition() + boxSize)));
		RayCollision collision;
		RayPlaneIntersection(CapsuleRay, Planes[i], collision);
		float p = (collision.collidedAt - ExtremeA).Length() - RadiusA;
		float distance = Maths::Clamp(p, 0.0f, CentreVector.Length());
		Vector3 ClosestPoint = ExtremeA + (CentreVector.Normalised() * distance);

		SphereVolume ClosestSphere = SphereVolume(RadiusA);
		Transform Transform = worldTransformA;
		Transform.SetPosition(ClosestPoint);

		if (AABBSphereIntersection(volumeB, worldTransformB, ClosestSphere, Transform, collisionInfo, true)) {
			collided = true;

			break;
		}
	}
	return collided;
}

bool CollisionDetection::OBBCapsuleIntersection(
	const CapsuleVolume& volumeA, const Transform& worldTransformA,
	const OBBVolume& volumeB, const Transform& worldTransformB, CollisionInfo& collisionInfo) {

	Matrix3 OBBTransform = Matrix3(worldTransformB.GetOrientation());
	Matrix3 invOBBTransform = Matrix3(worldTransformB.GetOrientation().Conjugate());

	Vector3 position = worldTransformB.GetPosition();
	Vector3 boxSize = volumeB.GetHalfDimensions();
	float RadiusA = volumeA.GetRadius();
	Vector3 ExtremeA = (worldTransformA.GetPosition() + worldTransformA.GetOrientation() * Vector3(0, 1, 0) * (volumeA.GetHalfHeight() - RadiusA));
	Vector3 ExtremeB = (worldTransformA.GetPosition() - worldTransformA.GetOrientation() * Vector3(0, 1, 0) * (volumeA.GetHalfHeight() - RadiusA));
	Vector3 CentreVector = ExtremeB - ExtremeA;

	Ray CapsuleRay((worldTransformA.GetPosition()) - CentreVector * 50, CentreVector);
	vector<Plane> Planes;
	vector<Vector3> normals = { OBBTransform * Vector3(1,0,0), OBBTransform * Vector3(0,1,0), OBBTransform * Vector3(0,0,1) };
	bool collided = false;

	for (int i = 0; i < normals.size(); i++) {
		Planes.push_back(Plane(normals[i], Vector3::Dot(-normals[i], worldTransformB.GetPosition() + boxSize)));
		RayCollision collision;
		RayPlaneIntersection(CapsuleRay, Planes[i], collision);
		float p = (collision.collidedAt - ExtremeA).Length() - RadiusA;
		float distance = Maths::Clamp(p, 0.0f, CentreVector.Length());
		Vector3 ClosestPoint = ExtremeA + (CentreVector.Normalised() * distance);

		SphereVolume ClosestSphere = SphereVolume(RadiusA);
		Transform Transform = worldTransformA;
		Transform.SetPosition(ClosestPoint);

		if (OBBSphereIntersection(volumeB, worldTransformB, ClosestSphere, Transform, collisionInfo, true)) {
			collided = true;

			break;
		}
	}
	return collided;
}

bool CollisionDetection::CapsuleCapsuleIntersection(const CapsuleVolume& volumeA, const Transform& worldTransformA,
	const CapsuleVolume& volumeB, const Transform& worldTransformB, CollisionInfo& collisionInfo) {
	Vector3 ExtremeA1 = (worldTransformA.GetPosition() + worldTransformA.GetOrientation() * Vector3(0, 1, 0) * (volumeA.GetHalfHeight() - volumeA.GetRadius()));
	Vector3 ExtremeB1 = (worldTransformA.GetPosition() - worldTransformA.GetOrientation() * Vector3(0, 1, 0) * (volumeA.GetHalfHeight() - volumeA.GetRadius()));
	Vector3 CentreVector1 = ExtremeB1 - ExtremeA1;

	Vector3 ExtremeA2 = (worldTransformB.GetPosition() + worldTransformB.GetOrientation() * Vector3(0, 1, 0) * (volumeB.GetHalfHeight() - volumeB.GetRadius()));
	Vector3 ExtremeB2 = (worldTransformB.GetPosition() - worldTransformB.GetOrientation() * Vector3(0, 1, 0) * (volumeB.GetHalfHeight() - volumeB.GetRadius()));
	Vector3 CentreVector2 = ExtremeB2 - ExtremeA2;

	Vector3 a = ExtremeA2 - ExtremeA1;
	Vector3 b = ExtremeB2 - ExtremeA1;
	Vector3 c = ExtremeA2 - ExtremeB1;
	Vector3 d = ExtremeB2 - ExtremeB1;

	float a2 = Vector3::Dot(a, a);
	float b2 = Vector3::Dot(b, b);
	float c2 = Vector3::Dot(c, c);
	float d2 = Vector3::Dot(d, d);

	Vector3 bestEnd1;
	if (c2 < a2 || c2 < b2 || d2 < a2 || d2 < b2) {
		bestEnd1 = ExtremeB1;
	}
	else {
		bestEnd1 = ExtremeA1;
	}

	Vector3 e = ExtremeA1 - ExtremeA2;
	Vector3 f = ExtremeB1 - ExtremeA2;
	Vector3 g = ExtremeA1 - ExtremeB2;
	Vector3 h = ExtremeB1 - ExtremeB2;

	float e2 = Vector3::Dot(a, a);
	float f2 = Vector3::Dot(b, b);
	float g2 = Vector3::Dot(c, c);
	float h2 = Vector3::Dot(d, d);

	Vector3 bestEnd2;
	if (g2 < e2 || g2 < f2 || h2 < e2 || h2 < f2) {
		bestEnd2 = ExtremeB2;
	}
	else {
		bestEnd2 = ExtremeA2;
	}

	float p = Vector3::Dot(bestEnd2 - ExtremeA1, CentreVector1) / Vector3::Dot(CentreVector1, CentreVector1);
	float distance = Maths::Clamp(p, 0.0f, 1.0f);
	Vector3 ClosestPointA = ExtremeA1 + (CentreVector1 * distance);

	float p2 = Vector3::Dot(bestEnd1 - ExtremeA2, CentreVector2) / Vector3::Dot(CentreVector2, CentreVector2);
	float distance2 = Maths::Clamp(p2, 0.0f, 1.0f);
	Vector3 ClosestPointB = ExtremeA2 + (CentreVector2 * distance2);

	SphereVolume ClosestSphereA = SphereVolume(volumeA.GetRadius());
	Transform TransformA = worldTransformA;
	TransformA.SetPosition(ClosestPointA);

	SphereVolume ClosestSphereB = SphereVolume(volumeB.GetRadius());
	Transform TransformB = worldTransformB;
	TransformB.SetPosition(ClosestPointB);

	return SphereIntersection(ClosestSphereA, TransformA, ClosestSphereB, TransformB, collisionInfo);
}