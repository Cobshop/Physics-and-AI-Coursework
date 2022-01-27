#pragma once
#include "Transform.h"
#include "CollisionVolume.h"
#include "AABBVolume.h"

#include "PhysicsObject.h"
#include "RenderObject.h"

#include <vector>

using std::vector;

namespace NCL {
	namespace CSC8503 {

		enum class ObjectType {
			Pull,
			Push,
			Spin,
			Crumble,
			Gravity,
			SurfaceChanger,
			Wall,
			Floor,
			Nothing,
		};

		class GameObject	{
		public:

			GameObject(string name = "");
			~GameObject();

			void SetBoundingVolume(CollisionVolume* vol) {
				boundingVolume = vol;
			}

			const CollisionVolume* GetBoundingVolume() const {
				return boundingVolume;
			}

			bool IsActive() const {
				return isActive;
			}

			bool IsTrigger() const {
				return isTrigger;
			}

			void SetTrigger(bool trigger) {
				isTrigger = trigger;
			}

			void SetType(ObjectType otype) {
				type = otype;
			}

			ObjectType GetType() {
				return type;
			}

			string GetTypeName() {
				switch (type) {
				case ObjectType::Pull: return "Pull"; break;
				case ObjectType::Push: return "Push"; break;
				case ObjectType::Spin: return "Spin"; break;
				case ObjectType::Crumble: return "Crumble"; break;
				case ObjectType::Gravity: return "Gravity"; break;
				case ObjectType::SurfaceChanger: return "SurfaceChanger"; break;
				case ObjectType::Wall: return "Wall"; break;
				case ObjectType::Floor: return "Floor"; break;
				case ObjectType::Nothing: return "Nothing"; break;
				}
			}

			void Deteriorate() {
				if (deteriorate) {
					renderObject->GetTransform()->SetScale(renderObject->GetTransform()->GetScale()
						- Vector3(renderObject->GetTransform()->GetScale().x/ 200, renderObject->GetTransform()->GetScale().y / 200, renderObject->GetTransform()->GetScale().z / 200));
					boundingVolume = (CollisionVolume*) new AABBVolume(renderObject->GetTransform()->GetScale()
						- Vector3(renderObject->GetTransform()->GetScale().x / 400, renderObject->GetTransform()->GetScale().y / 400, renderObject->GetTransform()->GetScale().z / 400));
				}
			}

			void ChangeSurface() {
				if (type == ObjectType::SurfaceChanger) {
					if(physicsObject->GetFriction() <= 0.9)
						physicsObject->SetFriction(physicsObject->GetFriction() + 0.1f);
					else
						physicsObject->SetFriction(physicsObject->GetFriction() - 0.9);
				}
			}

			Transform& GetTransform() {
				return transform;
			}

			RenderObject* GetRenderObject() const {
				return renderObject;
			}

			PhysicsObject* GetPhysicsObject() const {
				return physicsObject;
			}

			void SetRenderObject(RenderObject* newObject) {
				renderObject = newObject;

				initialScale = renderObject->GetTransform()->GetScale();
			}

			void SetPhysicsObject(PhysicsObject* newObject) {
				physicsObject = newObject;
			}

			const string& GetName() const {
				return name;
			}

			virtual void OnCollisionBegin(GameObject* otherObject) {
				//std::cout << "OnCollisionBegin event occured!\n";
			}

			virtual void OnCollisionEnd(GameObject* otherObject) {
				//std::cout << "OnCollisionEnd event occured!\n";
			}

			bool GetBroadphaseAABB(Vector3&outsize) const;

			void UpdateBroadphaseAABB();

			void SetWorldID(int newID) {
				worldID = newID;
			}

			int		GetWorldID() const {
				return worldID;
			}

			bool	deteriorate;
			float	surfaceType;
			Vector3 initialScale;
			bool	player;

		protected:
			Transform			transform;

			CollisionVolume*	boundingVolume;
			PhysicsObject*		physicsObject;
			RenderObject*		renderObject;

			bool	isActive;
			bool	isTrigger;
			
			int		worldID;
			string	name;
			ObjectType type;

			Vector3 broadphaseAABB;
		};
	}
}

