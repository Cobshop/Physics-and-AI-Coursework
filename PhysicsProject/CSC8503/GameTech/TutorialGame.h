#pragma once
#include "GameTechRenderer.h"
#include "../CSC8503Common/PhysicsSystem.h"
#include "../CSC8503Common/BehaviourAction.h"
#include "../CSC8503Common/BehaviourSequence.h"
#include "../CSC8503Common/BehaviourSelector.h"
#include "../CSC8503Common/NavigationGrid.h"
#include "StateGameObject.h"

namespace NCL {
	namespace CSC8503 {
		class TutorialGame		{
		public:
			TutorialGame(bool Game);
			~TutorialGame();

			virtual void UpdateGame(float dt);

			void Paused(float dt);
			void Menu(float dt);

			bool push;
			bool pushB;
			bool pop;
			bool GameEnd = false;

		protected:
			void InitialiseAssets();

			void InitCamera();
			void UpdateKeys();

			void InitWorld();

			void InitGameA();
			void InitGameB();
	
			bool SelectObject();
			void MoveSelectedObject();
			void DebugObjectMovement();
			void LockedObjectMovement();
			
			void CollectPowerUp();
			void CollectKey();
			void OpenChest();
			void UseHookshot();
			void ReactingBlocks();

			vector<Vector3> Pathfinding(GameObject* target);
			void RedBehaviourTree();
			Vector3 RedPathfinding(GameObject* target);
			void RedMove(GameObject* target);

			void GameAUpdates(float dt);
			void GameBUpdates(float dt);
			void DrawMaze(const std::string& filename);
			
			void GiveObjectInfo(GameObject* Object);

			GameObject* AddSphereToWorld(const Vector3& position, Vector4 colour, float radius, bool hollow, float inverseMass = 10.0f, float elasticity = 0.66, bool player = false);
			GameObject* AddCubeToWorld(const Vector3& position, Vector3 dimensions, Vector4 colour, float inverseMass = 10.0f, float elasticity = 0.66, bool OBB = false, float friction = 0.3f);
			GameObject* AddKeyToWorld(const Vector3& position, Vector3 dimensions, Vector4 colour, float inverseMass, float elasticity, float friction);
			GameObject* AddChestToWorld(const Vector3& position, Vector3 dimensions, Vector4 colour, float inverseMass, float elasticity, float friction, float angle);
			GameObject* AddCapsuleToWorld(const Vector3& position, float halfHeight, float radius, float inverseMass = 10.0f, float elasticity = 0.66, Vector4 colour = Vector4(1,1,1,1));
			GameObject* AddBonusToWorld(const Vector3& position, float radius);

			GameObject* TestPaths(vector<GameObject*> targets);

			void SetType(GameObject* object);

			StateGameObject* AddStateObjectToWorld(const Vector3& position, float radius);

			GameTechRenderer*	renderer;
			PhysicsSystem*		physics;
			GameWorld*			world;

			bool Controls;
			bool useGravity;
			bool inSelectionMode;
			bool Hookshot;
			bool usingHook;
			bool hookActive;
			bool changeGravity;
			bool Text = true;
			bool GameB = false;
			bool freeCam = false;
			int hookTimer;
			int confettiTimer;
			int hookDirection;
			float platformMover;
			float cameraMover;
			int PlayerPoints;
			int PlayerKeys;
			int EnemyPoints;
			int EnemyKeys;
			int BallType;

			std::chrono::high_resolution_clock::time_point tBegin;
			string Time;
			float TimeSec;

			Vector4 Ocolour;
			Vector3 RedAim;

			float		forceMagnitude;

			GameObject* selectionObject = nullptr;
			GameObject* zObject = nullptr;
			GameObject* Ball = nullptr;
			GameObject* PowerUp = nullptr;
			GameObject* Hook = nullptr;
			GameObject*	Confetti1 = nullptr;
			GameObject* Confetti2 = nullptr;
			PositionConstraint* hookConstraint = nullptr;

			BehaviourState state = Ongoing;
			BehaviourSequence* sequence = new BehaviourSequence("Room Sequence");
			BehaviourSelector* selection = new BehaviourSelector("Loot Selection");

			vector<GameObject*> movingPlatforms;
			vector<GameObject*> Keys;
			vector<GameObject*> Chests;
			GameObject* target = nullptr;

			GameObject* Red = nullptr;

			OGLMesh*	capsuleMesh = nullptr;
			OGLMesh*	cubeMesh	= nullptr;
			OGLMesh*	sphereMesh	= nullptr;
			OGLTexture* basicTex	= nullptr;
			OGLTexture* sphereTex = nullptr;
			OGLTexture* chestTex = nullptr;
			OGLShader*	basicShader = nullptr;

			//Coursework Meshes
			OGLMesh*	charMeshA	= nullptr;
			OGLMesh*	charMeshB	= nullptr;
			OGLMesh*	enemyMesh	= nullptr;
			OGLMesh*	bonusMesh	= nullptr;
			OGLMesh*	keyMesh = nullptr;
			OGLMesh*	chestMesh = nullptr;

			//Coursework Additional functionality	
			GameObject* lockedObject	= nullptr;
			Vector3 lockedOffset		= Vector3(0, 14, 20);
			void LockCameraToObject(GameObject* o) {
				lockedObject = o;
			}

		};
	}
}

