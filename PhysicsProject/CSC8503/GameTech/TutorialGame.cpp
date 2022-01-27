#include "TutorialGame.h"
#include "../CSC8503Common/GameWorld.h"
#include "../../Plugins/OpenGLRendering/OGLMesh.h"
#include "../../Plugins/OpenGLRendering/OGLShader.h"
#include "../../Plugins/OpenGLRendering/OGLTexture.h"
#include "../../Common/TextureLoader.h"
#include "../../Common/maths.h"
#include "../../Common/Assets.h"
#include "iostream"
#include <chrono>

#include <fstream>

#define USE_MATH_DEFINES

using namespace NCL;
using namespace CSC8503;

TutorialGame::TutorialGame(bool Game)	{
	world		= new GameWorld();
	renderer	= new GameTechRenderer(*world);
	physics		= new PhysicsSystem(*world);

	push = false;
	pushB = false;
	pop = false;
	forceMagnitude	= 10.0f;
	useGravity		= false;
	inSelectionMode = false;
	Hookshot = false;
	usingHook = false;
	hookActive = false;
	changeGravity = false;
	hookDirection = 1;
	hookTimer = 0;
	platformMover = 0;
	cameraMover = 0;
	BallType = 0;
	PlayerPoints = 0;
	EnemyPoints = 0;
	PlayerKeys = 0;
	EnemyKeys = 0;
	confettiTimer = 0;
	RedAim = Vector3(0, 0, 0);
	GameB = Game;

	Debug::SetRenderer(renderer);

	InitialiseAssets();
}

/*

Each of the little demo scenarios used in the game uses the same 2 meshes, 
and the same texture and shader. There's no need to ever load in anything else
for this module, even in the coursework, but you can add it if you like!

*/
void TutorialGame::InitialiseAssets() {
	auto loadFunc = [](const string& name, OGLMesh** into) {
		*into = new OGLMesh(name);
		(*into)->SetPrimitiveType(GeometryPrimitive::Triangles);
		(*into)->UploadToGPU();
	};

	loadFunc("cube.msh"		 , &cubeMesh);
	loadFunc("sphere.msh"	 , &sphereMesh);
	loadFunc("Male1.msh"	 , &charMeshA);
	loadFunc("courier.msh"	 , &charMeshB);
	loadFunc("security.msh"	 , &enemyMesh);
	loadFunc("coin.msh"		 , &bonusMesh);
	loadFunc("capsule.msh"	 , &capsuleMesh);
	loadFunc("key.msh"		 , &keyMesh);
	loadFunc("chest.msh", &chestMesh);

	basicTex	= (OGLTexture*)TextureLoader::LoadAPITexture("checkerboard_1.png");
	sphereTex = (OGLTexture*)TextureLoader::LoadAPITexture("checkerboard.png");
	chestTex = (OGLTexture*)TextureLoader::LoadAPITexture("chest.png");

	basicShader = new OGLShader("GameTechVert.glsl", "GameTechFrag.glsl");

	InitCamera();
	InitWorld();
}

TutorialGame::~TutorialGame()	{
	delete cubeMesh;
	delete sphereMesh;
	delete charMeshA;
	delete charMeshB;
	delete enemyMesh;
	delete bonusMesh;
	delete keyMesh;
	delete chestMesh;

	delete basicTex;
	delete basicShader;

	delete physics;
	delete renderer;
	delete world;
}

void TutorialGame::Paused(float dt) {
	
		Window::GetWindow()->ShowOSPointer(true);
		Window::GetWindow()->LockMouseToWindow(false);
		Vector2 screenSize = Window::GetWindow()->GetScreenSize();
		Vector2 screenMouse = Vector2(Window::GetMouse()->GetAbsolutePosition().x, Window::GetMouse()->GetAbsolutePosition().y);
		if (!GameEnd) {
			if (!Controls) {
				renderer->DrawString("Paused", Vector2(10, 20), Vector4(1, 1, 1, 1), 100);

				if (screenMouse.x > (125.0 / 1264.0)* screenSize.x && screenMouse.x < (407.0 / 1264.0)* screenSize.x && screenMouse.y > (251.0 / 681.0)* screenSize.y && screenMouse.y < (271.0 / 681.0)* screenSize.y) {
					renderer->DrawString("Controls", Vector2(10, 40), Vector4(1, 1, 0, 1), 30);

					if (Window::GetMouse()->ButtonPressed(NCL::MouseButtons::LEFT)) {
						Controls = true;
					}
				}
				else
					renderer->DrawString("Controls", Vector2(10, 40), Vector4(1, 1, 1, 1), 30);

				if (screenMouse.x > (126.0 / 1264.0)* screenSize.x && screenMouse.x < (551.0 / 1264.0)* screenSize.x && screenMouse.y > (387.0 / 681.0)* screenSize.y && screenMouse.y < (406.0 / 681.0)* screenSize.y) {
					renderer->DrawString("Exit to Menu", Vector2(10, 60), Vector4(1, 1, 0, 1), 30);

					if (Window::GetMouse()->ButtonPressed(NCL::MouseButtons::LEFT)) {
						pop = true;
					}
				}
				else
					renderer->DrawString("Exit to Menu", Vector2(10, 60), Vector4(1, 1, 1, 1), 30);

				if (screenMouse.x > (125.0 / 1264.0)* screenSize.x && screenMouse.x < (657.0 / 1264.0)* screenSize.x && screenMouse.y > (524.0 / 681.0)* screenSize.y && screenMouse.y < (546.0 / 681.0)* screenSize.y) {
					renderer->DrawString("Quit to Desktop", Vector2(10, 80), Vector4(1, 1, 0, 1), 30);

					if (Window::GetMouse()->ButtonPressed(NCL::MouseButtons::LEFT)) {
						push = true;
					}
				}
				else
					renderer->DrawString("Quit to Desktop", Vector2(10, 80), Vector4(1, 1, 1, 1), 30);
			}
			else {
				renderer->DrawString("O: Activate Hookshot", Vector2(10, 10));
				renderer->DrawString("U: Switch Hookshot Direction", Vector2(10, 20));
				renderer->DrawString("Space (Hold): Use Hookshot", Vector2(10, 30));
				renderer->DrawString("Left Click: Interact with World", Vector2(10, 40));

				if (screenMouse.x > (125.0 / 1264.0) * screenSize.x && screenMouse.x < (265.0 / 1264.0) * screenSize.x && screenMouse.y >(524.0 / 681.0) * screenSize.y && screenMouse.y < (546.0 / 681.0) * screenSize.y) {
					renderer->DrawString("Back", Vector2(10, 80), Vector4(1, 1, 0, 1), 30);

					if (Window::GetMouse()->ButtonPressed(NCL::MouseButtons::LEFT)) {
						Controls = false;
					}
				}
				else
					renderer->DrawString("Back", Vector2(10, 80), Vector4(1, 1, 1, 1), 30);
			}
		}
		else {
			if (GameB) {
				if (PlayerPoints > EnemyPoints) {
					renderer->DrawString("You Win!", Vector2(30, 20), Vector4(1, 1, 1, 1), 50);
					renderer->DrawString("Well Done!", Vector2(27, 40), Vector4(1, 1, 1, 1), 50);
					renderer->DrawString("Your Score was: " + std::to_string(PlayerPoints), Vector2(10, 60), Vector4(1, 1, 1, 1), 50);
				}
				else if (PlayerPoints < EnemyPoints) {
					renderer->DrawString("You Lost!", Vector2(30, 20), Vector4(1, 1, 1, 1), 50);
					renderer->DrawString("Well Done!", Vector2(27, 40), Vector4(1, 1, 1, 1), 50);
					renderer->DrawString("Your Score was: " + std::to_string(PlayerPoints), Vector2(10, 60), Vector4(1, 1, 1, 1), 50);
				}
				else if (PlayerPoints == EnemyPoints) {
					renderer->DrawString("You Draw!", Vector2(30, 20), Vector4(1, 1, 1, 1), 50);
					renderer->DrawString("Meh!", Vector2(40, 40), Vector4(1, 1, 1, 1), 50);
					renderer->DrawString("Your Score was: " + std::to_string(PlayerPoints), Vector2(10, 60), Vector4(1, 1, 1, 1), 50);
				}
			}
			else {
				renderer->DrawString("Congratulations!", Vector2(15, 20), Vector4(1, 1, 1, 1), 50);
				renderer->DrawString("Your Time Was " + Time, Vector2(5, 40), Vector4(1, 1, 1, 1), 50);

				if (TimeSec < 30) {
					renderer->DrawString("You Got a ", Vector2(10, 60), Vector4(1, 1, 1, 1), 50);
					renderer->DrawString("R", Vector2(60, 60), Vector4(1, 0, 0, 1), 50);
					renderer->DrawString("a", Vector2(65, 60), Vector4(1, 0.85, 0, 1), 50);
					renderer->DrawString("i", Vector2(70, 60), Vector4(1, 1, 0, 1), 50);
					renderer->DrawString("n", Vector2(75, 60), Vector4(0, 1, 0, 1), 50);
					renderer->DrawString("b", Vector2(80, 60), Vector4(0.5, 0.8, 0.9, 1), 50);
					renderer->DrawString("o", Vector2(85, 60), Vector4(0, 0, 1, 1), 50);
					renderer->DrawString("w", Vector2(90, 60), Vector4(1, 0, 1, 1), 50);
					renderer->DrawString("Time", Vector2(30, 80), Vector4(1, 1, 1, 1), 50);
				}
				else if (TimeSec < 60) {
					renderer->DrawString("You Got a ", Vector2(10, 60), Vector4(1, 1, 1, 1), 50);
					renderer->DrawString("Gold", Vector2(60, 60), Vector4(0.8, 0.7, 0.2, 1), 50);
					renderer->DrawString("Time", Vector2(40, 80), Vector4(1, 1, 1, 1), 50);
				}
				else if (TimeSec < 120) {
					renderer->DrawString("You Got a ", Vector2(10, 60), Vector4(1, 1, 1, 1), 50);
					renderer->DrawString("Silver", Vector2(60, 60), Vector4(0.5, 0.5, 0.5, 1), 50);
					renderer->DrawString("Time", Vector2(40, 80), Vector4(1, 1, 1, 1), 50);
				}
				else if (TimeSec < 300) {
					renderer->DrawString("You Got a ", Vector2(10, 60), Vector4(1, 1, 1, 1), 50);
					renderer->DrawString("Bronze", Vector2(60, 60), Vector4(0.8, 0.5, 0.2, 1), 50);
					renderer->DrawString("Time", Vector2(40, 80), Vector4(1, 1, 1, 1), 50);
				}
				else {
					renderer->DrawString("You Got a ", Vector2(10, 60), Vector4(1, 1, 1, 1), 50);
					renderer->DrawString("Bad", Vector2(60, 60), Vector4(0, 0, 0, 1), 50);
					renderer->DrawString("Time", Vector2(40, 80), Vector4(1, 1, 1, 1), 50);
				}
					
			}
			if (screenMouse.x > (125.0 / 1264.0)* screenSize.x && screenMouse.x < (265.0 / 1264.0)* screenSize.x && screenMouse.y > (524.0 / 681.0)* screenSize.y && screenMouse.y < (546.0 / 681.0)* screenSize.y) {
				renderer->DrawString("Back", Vector2(10, 80), Vector4(1, 1, 0, 1), 30);

				if (Window::GetMouse()->ButtonPressed(NCL::MouseButtons::LEFT)) {
					pop = true;
					GameEnd = false;
				}
			}
			else
				renderer->DrawString("Back", Vector2(10, 80), Vector4(1, 1, 1, 1), 30);
		}

	renderer->Update(dt);
	Debug::FlushRenderables(dt);
	renderer->Render();
}

void TutorialGame::Menu(float dt) {
	renderer->DrawString("GAME", Vector2(10, 20), Vector4(1,1,1,1), 100);

	Window::GetWindow()->ShowOSPointer(true);
	Window::GetWindow()->LockMouseToWindow(false);
	Vector2 screenSize = Window::GetWindow()->GetScreenSize();
	Vector2 screenMouse = Vector2(Window::GetMouse()->GetAbsolutePosition().x, Window::GetMouse()->GetAbsolutePosition().y);
	
	if (screenMouse.x > (126.0/1264.0) * screenSize.x && screenMouse.x < (515.0 / 1264.0) * screenSize.x && screenMouse.y > (251.0 / 681.0)* screenSize.y && screenMouse.y < (271.0 / 681.0)* screenSize.y) {
		renderer->DrawString("Play Game A", Vector2(10, 40), Vector4(1,1,0,1), 30);

		if (Window::GetMouse()->ButtonPressed(NCL::MouseButtons::LEFT)) {
			push = true;
		}
	}else
		renderer->DrawString("Play Game A", Vector2(10, 40), Vector4(1, 1, 1, 1), 30);

	if (screenMouse.x > (124.0 / 1264.0)* screenSize.x && screenMouse.x < (516.0 / 1264.0)* screenSize.x && screenMouse.y > (385.0 / 681.0)* screenSize.y && screenMouse.y < (407.0 / 681.0)* screenSize.y) {
		renderer->DrawString("Play Game B", Vector2(10, 60), Vector4(1, 1, 0, 1), 30);

		if (Window::GetMouse()->ButtonPressed(NCL::MouseButtons::LEFT)) {
			pushB = true;
		}
	}
	else
		renderer->DrawString("Play Game B", Vector2(10, 60), Vector4(1, 1, 1, 1), 30);

	if (screenMouse.x > (125.0 / 1264.0)* screenSize.x && screenMouse.x < (445.0 / 1264.0)* screenSize.x && screenMouse.y > (524.0 / 681.0)* screenSize.y && screenMouse.y < (544.0 / 681.0)* screenSize.y) {
		renderer->DrawString("Exit Game", Vector2(10, 80), Vector4(1, 1, 0, 1), 30);

		if (Window::GetMouse()->ButtonPressed(NCL::MouseButtons::LEFT)) {
			pop = true;
		}
	}
	else
		renderer->DrawString("Exit Game", Vector2(10, 80), Vector4(1, 1, 1, 1), 30);

	cameraMover += 0.05 * dt;
	
	world->GetMainCamera()->SetPosition(Vector3(70.5f, 200 + 200 * sin(cameraMover - PI/2), -8));

	platformMover += 0.5 * dt;
	for (int i = 0; i < movingPlatforms.size(); i++) {
		movingPlatforms[i]->GetRenderObject()->GetTransform()->SetPosition(Vector3(movingPlatforms[i]->GetRenderObject()->GetTransform()->GetPosition().x
			, movingPlatforms[i]->GetRenderObject()->GetTransform()->GetPosition().y, (80 - (movingPlatforms[i]->GetRenderObject()->GetTransform()->GetScale().z) / 2) * sin(2 * i + platformMover)));
	}
	if (platformMover >= 2 * PI)
		platformMover -= 2 * PI;

	if (cameraMover >= 2 * PI)
		cameraMover -= 2 * PI;

	world->GetMainCamera()->SetPitch(7);
	world->GetMainCamera()->SetYaw(92);

	renderer->Update(dt);
	Debug::FlushRenderables(dt);
	renderer->Render();
}

void TutorialGame::UpdateGame(float dt) {
	if (!inSelectionMode || freeCam) {
		world->GetMainCamera()->UpdateCamera(dt);
	}

	useGravity = true;
	physics->UseGravity(useGravity);

	SelectObject();
	MoveSelectedObject();
	physics->Update(dt);
	if (changeGravity)
		physics->SetGravity(Vector3(0, 0, 98.16));
	else
		physics->SetGravity(Vector3(0, -98.16, 0));

	if (lockedObject != nullptr) {
		Vector3 objPos = lockedObject->GetTransform().GetPosition();
		Vector3 camPos = objPos + lockedOffset;

		Matrix4 temp = Matrix4::BuildViewMatrix(camPos, objPos, Vector3(0,1,0));

		Matrix4 modelMat = temp.Inverse();

		Quaternion q(modelMat);
		Vector3 angles = q.ToEuler(); //nearly there now!

		world->GetMainCamera()->SetPosition(camPos);
		world->GetMainCamera()->SetPitch(angles.x);
		world->GetMainCamera()->SetYaw(angles.y);

		//Debug::DrawAxisLines(lockedObject->GetTransform().GetMatrix(), 2.0f);
	}

	if (Window::GetKeyboard()->KeyPressed(NCL::KeyboardKeys::M)) {
		freeCam = !freeCam;
	}

	world->UpdateWorld(dt);

	if (Ball && !GameB) {
		GameAUpdates(dt);
	}

	if (Ball && GameB) {
		GameBUpdates(dt);
	}
		
	
	renderer->Update(dt);
	

	Debug::FlushRenderables(dt);
	renderer->Render();
}

void TutorialGame::GameAUpdates(float dt) {
	inSelectionMode = true;

	platformMover += 0.5 * dt;

	for (int i = 0; i < movingPlatforms.size(); i++) {
		movingPlatforms[i]->GetRenderObject()->GetTransform()->SetPosition(Vector3(movingPlatforms[i]->GetRenderObject()->GetTransform()->GetPosition().x
			, movingPlatforms[i]->GetRenderObject()->GetTransform()->GetPosition().y, (80 - (movingPlatforms[i]->GetRenderObject()->GetTransform()->GetScale().z)/2) * sin(2 * i + platformMover)));
	}

	if (platformMover >= 2 * PI)
		platformMover -= 2 * PI;

	Ball->GetTransform().SetPosition(Vector3(-80, Ball->GetTransform().GetPosition().y, Ball->GetTransform().GetPosition().z));

	if (!freeCam) {
		world->GetMainCamera()->SetPosition(Vector3(70.5f, Ball->GetTransform().GetPosition().y, -8));
		world->GetMainCamera()->SetPitch(7);
		world->GetMainCamera()->SetYaw(92);
	}

	CollectPowerUp();
	if (Hookshot)
		UseHookshot();

	confettiTimer++;

	if (Confetti1) {
		Confetti1->GetRenderObject()->GetTransform()->SetPosition(Vector3(-80, Confetti1->GetRenderObject()->GetTransform()->GetPosition().y, Confetti1->GetRenderObject()->GetTransform()->GetPosition().z));
	}
	if (Confetti2) {
		Confetti2->GetRenderObject()->GetTransform()->SetPosition(Vector3(-80, Confetti2->GetRenderObject()->GetTransform()->GetPosition().y, Confetti2->GetRenderObject()->GetTransform()->GetPosition().z));
	}

	if (confettiTimer > 400) {
		if (Confetti1)
			world->RemoveGameObject(Confetti1);
		if (Confetti2)
			world->RemoveGameObject(Confetti2);

		BallType++;
		if (BallType > 2) {
			BallType = 0;
		}
		SphereVolume* volume = new SphereVolume(5);
		CapsuleVolume* volume2 = new CapsuleVolume(5, 5);
		if (BallType == 0) {
			Confetti1 = AddSphereToWorld(Vector3(-80, 400, (30 + rand() % 30)), Vector4(1, 1, 1, 1), 5, false, 400.0f, 0.66f);
			Confetti2 = AddCapsuleToWorld(Vector3(-80, 400, -(30 + rand() % 30)), 8, 4, 400.0f);
			Confetti2->GetRenderObject()->GetTransform()->SetOrientation(Quaternion::AxisAngleToQuaterion(Vector3(1, 0, 0), 45));
		}
		if (BallType == 1) {
			Confetti1 = AddCubeToWorld(Vector3(-80, 400, -(30 + rand() % 30)), Vector4(5, 5, 5, 5), Vector4(1, 1, 1, 1), 400.0f, 0.2f, true);
			Confetti2 = AddCapsuleToWorld(Vector3(-80, 400, (30 + rand() % 30)), 8, 4, 400.0f);
			Confetti2->GetRenderObject()->GetTransform()->SetOrientation(Quaternion::AxisAngleToQuaterion(Vector3(1, 0, 0), 45));
		}
		if (BallType == 2) {
			Confetti1 = AddCubeToWorld(Vector3(-80, 400, -(30 + rand() % 30)), Vector4(5, 5, 5, 5), Vector4(1, 1, 1, 1), 400.0f, 0.2f, true);
			Confetti2 = AddSphereToWorld(Vector3(-80, 400, (30 + rand() % 30)), Vector4(1, 1, 1, 1), 5, false, 400.0f, 0.8f);
		}
		confettiTimer = 0;
	}

	if (pop == false && Ball->GetTransform().GetPosition().y > 360) {
		GameEnd = true;
	}
	auto tEnd = std::chrono::high_resolution_clock::now();
	Time = std::to_string(std::chrono::floor<std::chrono::minutes>(std::chrono::duration<float>(tEnd - tBegin)).count()) + "m "
		+ std::to_string(std::chrono::round<std::chrono::seconds>(std::chrono::duration<float>(tEnd - tBegin)).count() % 60) + "s ";
	TimeSec = std::chrono::round<std::chrono::seconds>(std::chrono::duration<float>(tEnd - tBegin)).count();
	renderer->DrawString("Time: " + Time, Vector2(60, 80), Vector4(1, 1, 1, 1), 25);
}

void TutorialGame::GameBUpdates(float dt) {
	if (Chests.size() > 0 && Keys.size() > 0) {
		if (state == Ongoing) {
			state = selection->Execute(1.0f); //fake dt
			RedMove(target);
		}
		if (state == Success) {
			selection->Reset();
			state = Initialise;
			std::cout << "Woo Points";
			RedAim = RedPathfinding(target);
			state = selection->Execute(1.0f);
		}
		if (state == Failure) {
			selection->Reset();
			state = Initialise;
			state = selection->Execute(1.0f);
		}
	}
	else if(pop == false){
		GameEnd = true;
	}	

	if (abs(Ball->GetPhysicsObject()->GetLinearVelocity().z) < 50) {
		if (Window::GetKeyboard()->KeyDown(KeyboardKeys::UP)) {
			Ball->GetPhysicsObject()->SetLinearVelocity(Vector3(Ball->GetPhysicsObject()->GetLinearVelocity().x + 1,
																Ball->GetPhysicsObject()->GetLinearVelocity().y,
																Ball->GetPhysicsObject()->GetLinearVelocity().z));
		}
		if (Window::GetKeyboard()->KeyDown(KeyboardKeys::DOWN)) {
			Ball->GetPhysicsObject()->SetLinearVelocity(Vector3(Ball->GetPhysicsObject()->GetLinearVelocity().x - 1,
																Ball->GetPhysicsObject()->GetLinearVelocity().y,
																Ball->GetPhysicsObject()->GetLinearVelocity().z));
		}
	}
	if (abs(Ball->GetPhysicsObject()->GetLinearVelocity().x) < 50) {
		if (Window::GetKeyboard()->KeyDown(KeyboardKeys::LEFT)) {
			Ball->GetPhysicsObject()->SetLinearVelocity(Vector3(Ball->GetPhysicsObject()->GetLinearVelocity().x,
																Ball->GetPhysicsObject()->GetLinearVelocity().y,
																Ball->GetPhysicsObject()->GetLinearVelocity().z - 1));
		}
		if (Window::GetKeyboard()->KeyDown(KeyboardKeys::RIGHT)) {
			Ball->GetPhysicsObject()->SetLinearVelocity(Vector3(Ball->GetPhysicsObject()->GetLinearVelocity().x,
																Ball->GetPhysicsObject()->GetLinearVelocity().y,
																Ball->GetPhysicsObject()->GetLinearVelocity().z + 1));
		}
	}

	renderer->DrawString("PlayerPoints: " + std::to_string(PlayerPoints), Vector2(10, 10), Vector4(1, 1, 1, 1), 25);
	renderer->DrawString("EnemyPoints: " + std::to_string(EnemyPoints), Vector2(60, 10), Vector4(1,0,0,1), 25);
	renderer->DrawString("PlayerKeys: " + std::to_string(PlayerKeys), Vector2(10, 20), Vector4(1, 1, 1, 1), 25);
	renderer->DrawString("EnemyKeys: " + std::to_string(EnemyKeys), Vector2(60, 20), Vector4(1, 0, 0, 1), 25);

	CollectKey();
	OpenChest();

	if (!freeCam) {
		world->GetMainCamera()->SetPosition(Vector3(Ball->GetTransform().GetPosition().x, Ball->GetTransform().GetPosition().y + 360, Ball->GetTransform().GetPosition().z));
		world->GetMainCamera()->SetPitch(-90);
		world->GetMainCamera()->SetYaw(270);
	}

	for (int i = 0; i < Keys.size(); i++) {
		CollisionDetection::CollisionInfo info;
		if (CollisionDetection::ObjectIntersection(Red, Keys[i], info)) {
			world->RemoveGameObject(Keys[i]);
			Keys.erase(std::remove(Keys.begin(), Keys.end(), Keys[i]), Keys.end());
			EnemyKeys++;
		}
	}
}

void TutorialGame::UpdateKeys() {
	if (Window::GetKeyboard()->KeyPressed(KeyboardKeys::F1)) {
		InitWorld(); //We can reset the simulation at any time with F1
		selectionObject = nullptr;
		lockedObject	= nullptr;
	}

	if (Window::GetKeyboard()->KeyPressed(KeyboardKeys::F2)) {
		InitCamera(); //F2 will reset the camera to a specific default place
	}

	if (Window::GetKeyboard()->KeyPressed(KeyboardKeys::G)) {
		useGravity = !useGravity; //Toggle gravity!
		physics->UseGravity(useGravity);
	}
	//Running certain physics updates in a consistent order might cause some
	//bias in the calculations - the same objects might keep 'winning' the constraint
	//allowing the other one to stretch too much etc. Shuffling the order so that it
	//is random every frame can help reduce such bias.
	if (Window::GetKeyboard()->KeyPressed(KeyboardKeys::F9)) {
		world->ShuffleConstraints(true);
	}
	if (Window::GetKeyboard()->KeyPressed(KeyboardKeys::F10)) {
		world->ShuffleConstraints(false);
	}

	if (Window::GetKeyboard()->KeyPressed(KeyboardKeys::F7)) {
		world->ShuffleObjects(true);
	}
	if (Window::GetKeyboard()->KeyPressed(KeyboardKeys::F8)) {
		world->ShuffleObjects(false);
	}

	if (lockedObject) {
		LockedObjectMovement();
	}
	else {
		DebugObjectMovement();
	}
}

void TutorialGame::LockedObjectMovement() {
	Matrix4 view		= world->GetMainCamera()->BuildViewMatrix();
	Matrix4 camWorld	= view.Inverse();

	Vector3 rightAxis = Vector3(camWorld.GetColumn(0)); //view is inverse of model!

	//forward is more tricky -  camera forward is 'into' the screen...
	//so we can take a guess, and use the cross of straight up, and
	//the right axis, to hopefully get a vector that's good enough!

	Vector3 fwdAxis = Vector3::Cross(Vector3(0, 1, 0), rightAxis);
	fwdAxis.y = 0.0f;
	fwdAxis.Normalise();

	Vector3 charForward  = lockedObject->GetTransform().GetOrientation() * Vector3(0, 0, 1);
	Vector3 charForward2 = lockedObject->GetTransform().GetOrientation() * Vector3(0, 0, 1);

	float force = 100.0f;

	if (Window::GetKeyboard()->KeyDown(KeyboardKeys::LEFT)) {
		lockedObject->GetPhysicsObject()->AddForce(-rightAxis * force);
	}

	if (Window::GetKeyboard()->KeyDown(KeyboardKeys::RIGHT)) {
		Vector3 worldPos = selectionObject->GetTransform().GetPosition();
		lockedObject->GetPhysicsObject()->AddForce(rightAxis * force);
	}

	if (Window::GetKeyboard()->KeyDown(KeyboardKeys::UP)) {
		lockedObject->GetPhysicsObject()->AddForce(fwdAxis * force);
	}

	if (Window::GetKeyboard()->KeyDown(KeyboardKeys::DOWN)) {
		lockedObject->GetPhysicsObject()->AddForce(-fwdAxis * force);
	}

	if (Window::GetKeyboard()->KeyDown(KeyboardKeys::NEXT)) {
		lockedObject->GetPhysicsObject()->AddForce(Vector3(0,-10,0));
	}
}

void TutorialGame::DebugObjectMovement() {
//If we've selected an object, we can manipulate it with some key presses
	if (inSelectionMode && selectionObject) {
		//Twist the selected object!
		if (Window::GetKeyboard()->KeyDown(KeyboardKeys::LEFT)) {
			selectionObject->GetPhysicsObject()->AddTorque(Vector3(-10, 0, 0));
		}

		if (Window::GetKeyboard()->KeyDown(KeyboardKeys::RIGHT)) {
			selectionObject->GetPhysicsObject()->AddTorque(Vector3(10, 0, 0));
		}

		if (Window::GetKeyboard()->KeyDown(KeyboardKeys::NUM7)) {
			selectionObject->GetPhysicsObject()->AddTorque(Vector3(0, 10, 0));
		}

		if (Window::GetKeyboard()->KeyDown(KeyboardKeys::NUM8)) {
			selectionObject->GetPhysicsObject()->AddTorque(Vector3(0, -10, 0));
		}

		if (Window::GetKeyboard()->KeyDown(KeyboardKeys::RIGHT)) {
			selectionObject->GetPhysicsObject()->AddTorque(Vector3(10, 0, 0));
		}

		if (Window::GetKeyboard()->KeyDown(KeyboardKeys::UP)) {
			selectionObject->GetPhysicsObject()->AddForce(Vector3(0, 0, -10));
		}

		if (Window::GetKeyboard()->KeyDown(KeyboardKeys::DOWN)) {
			selectionObject->GetPhysicsObject()->AddForce(Vector3(0, 0, 10));
		}

		if (Window::GetKeyboard()->KeyDown(KeyboardKeys::NUM5)) {
			selectionObject->GetPhysicsObject()->AddForce(Vector3(0, -10, 0));
		}
	}

}

void TutorialGame::InitCamera() {
	world->GetMainCamera()->SetNearPlane(0.1f);
	world->GetMainCamera()->SetFarPlane(500.0f);
	world->GetMainCamera()->SetPitch(-15.0f);
	world->GetMainCamera()->SetYaw(315.0f);
	world->GetMainCamera()->SetPosition(Vector3(-60, 40, 60));
	lockedObject = nullptr;
}

void TutorialGame::InitWorld() {
	world->ClearAndErase();
	physics->Clear();

	if(!GameB)
		InitGameA();
	else
		InitGameB();
}

/*

Builds a game object that uses a sphere mesh for its graphics, and a bounding sphere for its
rigid body representation. This and the cube function will let you build a lot of 'simple' 
physics worlds. You'll probably need another function for the creation of OBB cubes too.

*/
GameObject* TutorialGame::AddSphereToWorld(const Vector3& position, Vector4 colour, float radius, bool hollow, float inverseMass, float elasticity, bool player) {
	GameObject* sphere = new GameObject();

	Vector3 sphereSize = Vector3(radius, radius, radius);
	SphereVolume* volume = new SphereVolume(radius);
	sphere->SetBoundingVolume((CollisionVolume*)volume);

	sphere->GetTransform()
		.SetScale(sphereSize)
		.SetPosition(position);

	sphere->SetRenderObject(new RenderObject(&sphere->GetTransform(), sphereMesh, sphereTex, basicShader));
	sphere->SetPhysicsObject(new PhysicsObject(&sphere->GetTransform(), sphere->GetBoundingVolume()));
	sphere->player = player;

	sphere->GetRenderObject()->SetColour(colour);
	sphere->GetPhysicsObject()->SetInverseMass(inverseMass);
	sphere->GetPhysicsObject()->SetElasticity(elasticity);
	if(hollow)
		sphere->GetPhysicsObject()->InitHollowSphereInertia();
	else
		sphere->GetPhysicsObject()->InitSphereInertia();
	SetType(sphere);

	world->AddGameObject(sphere);

	return sphere;
}

GameObject* TutorialGame::AddCapsuleToWorld(const Vector3& position, float halfHeight, float radius, float inverseMass, float elasticity, Vector4 colour) {
	GameObject* capsule = new GameObject();

	CapsuleVolume* volume = new CapsuleVolume(halfHeight, radius);
	capsule->SetBoundingVolume((CollisionVolume*)volume);

	capsule->GetTransform()
		.SetScale(Vector3(radius* 2, halfHeight, radius * 2))
		.SetPosition(position);

	capsule->SetRenderObject(new RenderObject(&capsule->GetTransform(), capsuleMesh, basicTex, basicShader));
	capsule->SetPhysicsObject(new PhysicsObject(&capsule->GetTransform(), capsule->GetBoundingVolume()));

	capsule->GetPhysicsObject()->SetInverseMass(inverseMass);
	capsule->GetPhysicsObject()->SetElasticity(elasticity);
	capsule->GetPhysicsObject()->InitCubeInertia();
	capsule->GetRenderObject()->SetColour(colour);

	world->AddGameObject(capsule);

	return capsule;

}

GameObject* TutorialGame::AddCubeToWorld(const Vector3& position, Vector3 dimensions, Vector4 colour, float inverseMass, float elasticity, bool OBB, float friction) {
	GameObject* cube = new GameObject();

	if (OBB)
		cube->SetBoundingVolume((CollisionVolume*)new OBBVolume(dimensions));
	else
		cube->SetBoundingVolume((CollisionVolume*)new AABBVolume(dimensions));


	cube->GetTransform()
		.SetPosition(position)
		.SetScale(dimensions * 2);

	cube->SetRenderObject(new RenderObject(&cube->GetTransform(), cubeMesh, basicTex, basicShader));
	cube->SetPhysicsObject(new PhysicsObject(&cube->GetTransform(), cube->GetBoundingVolume()));

	cube->GetPhysicsObject()->SetInverseMass(inverseMass);
	cube->GetPhysicsObject()->SetElasticity(elasticity);
	cube->GetPhysicsObject()->SetFriction(friction);
	cube->GetPhysicsObject()->InitCubeInertia();
	cube->GetRenderObject()->SetColour(colour);
	SetType(cube);

	world->AddGameObject(cube);

	return cube;
}

GameObject* TutorialGame::AddKeyToWorld(const Vector3& position, Vector3 dimensions, Vector4 colour, float inverseMass, float elasticity, float friction) {
	GameObject* key = new GameObject();

	key->SetBoundingVolume((CollisionVolume*)new AABBVolume(dimensions * 10));

	key->GetTransform()
		.SetPosition(position)
		.SetScale(dimensions * 2);

	key->SetRenderObject(new RenderObject(&key->GetTransform(), keyMesh, sphereTex, basicShader));
	key->SetPhysicsObject(new PhysicsObject(&key->GetTransform(), key->GetBoundingVolume()));

	key->GetPhysicsObject()->SetInverseMass(inverseMass);
	key->GetPhysicsObject()->SetElasticity(elasticity);
	key->GetPhysicsObject()->SetFriction(friction);
	key->GetPhysicsObject()->InitCubeInertia();
	key->GetRenderObject()->SetColour(colour);
	key->SetTrigger(true);

	world->AddGameObject(key);

	return key;

}

GameObject* TutorialGame::AddChestToWorld(const Vector3& position, Vector3 dimensions, Vector4 colour, float inverseMass, float elasticity, float friction, float angle) {
	GameObject* key = new GameObject();

	key->SetBoundingVolume((CollisionVolume*)new AABBVolume(dimensions));

	key->GetTransform()
		.SetPosition(position)
		.SetScale(dimensions * 2);

	key->SetRenderObject(new RenderObject(&key->GetTransform(), chestMesh, chestTex, basicShader));
	key->SetPhysicsObject(new PhysicsObject(&key->GetTransform(), key->GetBoundingVolume()));

	key->GetPhysicsObject()->SetInverseMass(inverseMass);
	key->GetPhysicsObject()->SetElasticity(elasticity);
	key->GetPhysicsObject()->SetFriction(friction);
	key->GetPhysicsObject()->InitCubeInertia();
	key->GetRenderObject()->SetColour(colour);
	key->GetRenderObject()->GetTransform()->SetOrientation(Quaternion::AxisAngleToQuaterion(Vector3(0, 1, 0), angle));

	world->AddGameObject(key);

	return key;

}

void TutorialGame::SetType(GameObject* object) {
	if (object->GetRenderObject()->GetColour() == Vector4(1, 0, 0, 1))
		object->SetType(ObjectType::Push);
	
	if (object->GetRenderObject()->GetColour() == Vector4(0, 1, 1, 1))
		object->SetType(ObjectType::Pull);

	if (object->GetRenderObject()->GetColour() == Vector4(0.5, 0.5, 0.5, 1))
		object->SetType(ObjectType::Crumble);

	if (object->GetRenderObject()->GetColour() == Vector4(1, 0.5, 0.5, 1))
		object->SetType(ObjectType::Spin);

	if (object->GetRenderObject()->GetColour() == Vector4(0.6, 0.3, 1, 1))
		object->SetType(ObjectType::Gravity);

	if (object->GetRenderObject()->GetColour() == Vector4(1, 0.2, 0, 1))
		object->SetType(ObjectType::SurfaceChanger);

	if (object->GetRenderObject()->GetColour() == Vector4(0, 0, 1, 1))
		object->SetType(ObjectType::Floor);

	if (object->GetRenderObject()->GetColour() == Vector4(1, 1, 0, 1))
		object->SetType(ObjectType::Wall);

	if (object->GetRenderObject()->GetColour() == Vector4(1, 1, 1, 1))
		object->SetType(ObjectType::Nothing);
}

void TutorialGame::InitGameA() {
	tBegin = std::chrono::high_resolution_clock::now();

	AddCubeToWorld(Vector3(-100, 180, 0), Vector3(10, 180, 100), Vector4(1, 0, 1, 1), 0);
	for (int i = 0; i < 5; i++) {
		AddCubeToWorld(Vector3(-80, 50 + 20* i, -90), Vector3(10, 10, 10), Vector4(1, 1, 0, 1), 0);
	}
	for (int i = 0; i < 10; i++) {
		AddCubeToWorld(Vector3(-80, 170 + 20 * i, -90), Vector3(10, 10, 10), Vector4(1, 1, 0, 1), 0);
	}
	AddCubeToWorld(Vector3(-80, 10, -90), Vector3(10, 10, 10), Vector4(1, 1, 0, 1), 0);
	for (int i = 0; i < 18; i++) {
		AddCubeToWorld(Vector3(-80, 10 + 20 * i, 90), Vector3(10, 10, 10), Vector4(1, 1, 0, 1), 0);
	}

	AddCubeToWorld(Vector3(-80, 30, -90), Vector3(10, 10, 10), Vector4(0, 1, 1, 1), 0);
	AddCubeToWorld(Vector3(-80, 265, -30), Vector3(10, 10, 10), Vector4(0, 1, 1, 1), 0);

	AddCubeToWorld(Vector3(-80, 10, 30), Vector3(10, 10, 50), Vector4(0, 0, 1, 1), 0);
	AddCubeToWorld(Vector3(-80, 130, -30), Vector3(10, 10, 30), Vector4(0, 0, 1, 1), 0);
	AddCubeToWorld(Vector3(-80, 150, 85), Vector3(9, 50, 10), Vector4(1, 0.2, 0, 1), 0, 0.66f, true)->GetRenderObject()->GetTransform()->SetOrientation(Quaternion::AxisAngleToQuaterion(Vector3(1, 0, 0), 5));
	AddCubeToWorld(Vector3(-80, 300, 70), Vector3(9, 30, 10), Vector4(1, 0.2, 0, 1), 0, 0.66f, true)->GetRenderObject()->GetTransform()->SetOrientation(Quaternion::AxisAngleToQuaterion(Vector3(1, 0, 0), 45));
	AddCubeToWorld(Vector3(-80, 200, -10), Vector3(10, 10, 30), Vector4(0.5, 0.5, 0.5, 1), 0);

	AddCubeToWorld(Vector3(-80, 150, -90), Vector3(10, 10, 10), Vector4(0, 1, 1, 1), 0);
	AddCubeToWorld(Vector3(-80, 130, -70), Vector3(10, 10, 10), Vector4(1, 0, 0, 1), 0, 1);

	AddCubeToWorld(Vector3(-80, 10, -30), Vector3(10, 10, 10), Vector4(1, 0, 0, 1), 0, 1);
	AddCubeToWorld(Vector3(-80, 10, -50), Vector3(10, 10, 10), Vector4(1, 0, 0, 1), 0, 1);
	AddCubeToWorld(Vector3(-80, 10, -70), Vector3(10, 10, 10), Vector4(1, 0, 0, 1), 0, 1);

	movingPlatforms.push_back(AddCubeToWorld(Vector3(-80, 250, 0), Vector3(10, 5, 10), Vector4(1, 0.5, 0.5, 1), 0, 0.66f, true));
	//movingPlatforms.push_back(AddCubeToWorld(Vector3(-80, 350, 0), Vector3(10, 5, 30), Vector4(1, 0.5, 0.5, 1), 0, 0.66f, true));
	movingPlatforms.push_back(AddCapsuleToWorld(Vector3(-80, 350, 0), 10, 5, 0, 0.66f));
	movingPlatforms[1]->GetRenderObject()->GetTransform()->SetOrientation(Quaternion::AxisAngleToQuaterion(Vector3(1, 0, 0), 90));

	AddCubeToWorld(Vector3(-80, 70, 0), Vector3(10, 5, 30), Vector4(1, 0.5, 0.5, 1), 0, 0.66f, true);

	AddSphereToWorld(Vector3(-80, 220, 50), Vector4(0.6, 0.3, 1, 1), 5, false, 0, 0.66f);
	AddCubeToWorld(Vector3(-80, 200, 30), Vector3(10, 10, 10), Vector4(1, 0, 0, 1), 0, 1);

	Ball = AddSphereToWorld(Vector3(-80, 30, 70), Vector4(1, 1, 1, 1), 5, false, 400.0f, 0.66f, true);
	PowerUp = AddBonusToWorld(Vector3(-80, 100, -70), 1);
}

void TutorialGame::InitGameB() {
	AddCubeToWorld(Vector3(200, 0, 200), Vector3(210, 10, 210), Vector4(0, 1, 0, 1), 0);

	try {
		DrawMaze("TestGrid1.txt");
	}
	catch (const std::invalid_argument& iae) {
		std::cout << "unable to read data: " << iae.what() << "\n";
		exit(1);
	}

	Ball = AddSphereToWorld(Vector3(20, 20, 20), Vector4(1, 1, 1, 1), 5, false, 600.0f, 0.66f, true);//AddCubeToWorld(Vector3(-170, 20, 170), Vector4(5, 5, 5, 5), Vector4(1, 1, 1, 1), 400.0f, 0.2f, true);

	Red = AddSphereToWorld(Vector3(380, 20, 380), Vector4(1, 0, 0, 1), 5, false, 600.0f, 0.66f, true);

	Chests.push_back(AddChestToWorld(Vector3(20, 10, 145), Vector3(10, 10, 10), Vector4(1, 1, 1, 1), 0, 0.66f, 0.3f, 180));
	Chests.push_back(AddChestToWorld(Vector3(380, 10, 60), Vector3(10, 10, 10), Vector4(1, 1, 1, 1), 0, 0.66f, 0.3f, 180));
	Chests.push_back(AddChestToWorld(Vector3(260, 10, 100), Vector3(10, 10, 10), Vector4(1, 1, 1, 1), 0, 0.66f, 0.3f, 180));
	Chests.push_back(AddChestToWorld(Vector3(180, 10, 260), Vector3(10, 10, 10), Vector4(1, 1, 1, 1), 0, 0.66f, 0.3f, 90));
	Chests.push_back(AddChestToWorld(Vector3(300, 10, 380), Vector3(10, 10, 10), Vector4(1, 1, 1, 1), 0, 0.66f, 0.3f, 180));
	Chests.push_back(AddChestToWorld(Vector3(340, 10, 300), Vector3(10, 10, 10), Vector4(1, 1, 1, 1), 0, 0.66f, 0.3f, 0));

	RedAim = RedPathfinding(Keys[0]);

	RedBehaviourTree();
}

void TutorialGame::DrawMaze(const std::string& filename) {
	int nodeSize;
	int	gridWidth;
	int	gridHeight;

	std::ifstream infile(Assets::DATADIR + filename);

	string temp;
	int j = 0;

	if (infile.fail())
		throw std::invalid_argument("no file exists " + filename);

	for (temp; getline(infile, temp);) {

		if (temp[0] == 'x') {
			for (int i = 0; i < temp.length(); i++) {
				int k = rand() % 20;
				if(temp[i] == 'x')
					AddCubeToWorld(Vector3(i * 20, 20, j * 20), Vector3(10, 10, 10), Vector4(0, 0, 1, 1), 0);
				else if(k == 2)
					Keys.push_back(AddKeyToWorld(Vector3(i * 20 + 5, 15, j * 20 + 5), Vector3(1, 1, 1), Vector4(1, 1, 0, 1), 0, 0.66f, 0.3f));
			}
			j++;
		}
	}

	infile.close();
}

GameObject* TutorialGame::AddBonusToWorld(const Vector3& position, float radius) {
	GameObject* apple = new GameObject();

	SphereVolume* volume = new SphereVolume(4.8 * radius);
	apple->SetBoundingVolume((CollisionVolume*)volume);
	apple->GetTransform()
		.SetScale(Vector3(radius, radius, radius))
		.SetPosition(position);

	apple->GetTransform().SetOrientation(Quaternion(0, sin(PI / 4), 0, cos(PI / 4)));

	apple->SetRenderObject(new RenderObject(&apple->GetTransform(), bonusMesh, nullptr, basicShader));
	apple->SetPhysicsObject(new PhysicsObject(&apple->GetTransform(), apple->GetBoundingVolume()));

	apple->GetPhysicsObject()->SetInverseMass(0.0f);
	apple->SetTrigger(true);
	apple->GetPhysicsObject()->InitSphereInertia();

	world->AddGameObject(apple);

	return apple;
}

void TutorialGame::CollectPowerUp() {
	CollisionDetection::CollisionInfo info;
	if (CollisionDetection::ObjectIntersection(Ball, PowerUp, info)) {
		Hookshot = true;
		world->RemoveGameObject(PowerUp);
	}
}

void TutorialGame::CollectKey() {
	for (int i = 0; i < Keys.size(); i++) {
		CollisionDetection::CollisionInfo info;
		if (CollisionDetection::ObjectIntersection(Ball, Keys[i], info)) {
			world->RemoveGameObject(Keys[i]);
			Keys.erase(std::remove(Keys.begin(), Keys.end(), Keys[i]), Keys.end());
			PlayerKeys++;
		}
	}
}

void TutorialGame::OpenChest() {
	for (int i = 0; i < Chests.size(); i++) {
		CollisionDetection::CollisionInfo info;
		if (CollisionDetection::ObjectIntersection(Ball, Chests[i], info)) {
			if (PlayerKeys > 0) {
				world->RemoveGameObject(Chests[i]);
				Chests.erase(std::remove(Chests.begin(), Chests.end(), Chests[i]), Chests.end());
				PlayerKeys--;
				PlayerPoints+=15;
			}
		}
	}
}

StateGameObject* TutorialGame::AddStateObjectToWorld(const Vector3& position, float radius) {
	StateGameObject* apple = new StateGameObject();

	SphereVolume* volume = new SphereVolume(4.8 * radius);
	apple->SetBoundingVolume((CollisionVolume*)volume);
	apple->GetTransform()
		.SetScale(Vector3(radius, radius, radius))
		.SetPosition(position);

	apple->GetTransform().SetOrientation(Quaternion(0, sin(PI / 4), 0, cos(PI / 4)));

	apple->SetRenderObject(new RenderObject(&apple->GetTransform(), bonusMesh, nullptr, basicShader));
	apple->SetPhysicsObject(new PhysicsObject(&apple->GetTransform(), apple->GetBoundingVolume()));

	apple->GetPhysicsObject()->SetInverseMass(1.0f);
	apple->SetTrigger(true);
	apple->GetPhysicsObject()->InitSphereInertia();

	world->AddGameObject(apple);

	return apple;
}

void TutorialGame::UseHookshot() {
	Vector3 gravityChange = Vector3(1, 1, 1);;

	if (Window::GetKeyboard()->KeyPressed(NCL::KeyboardKeys::O)) {
		Text = false;
		usingHook = !usingHook;
	}
	
	if(Text) {
		renderer->DrawString("You found the hookshot! Press O to Activate", Vector2(5, 50));
	}

	if(Window::GetKeyboard()->KeyPressed(NCL::KeyboardKeys::U)) {
		hookDirection = -hookDirection;
	}
	
	if(changeGravity)
		gravityChange = Vector3(1, hookDirection, -1);
	else
		gravityChange = Vector3(1, 1, hookDirection);

	if (usingHook) {
		Debug::DrawLine(Ball->GetTransform().GetPosition(), Ball->GetTransform().GetPosition() + gravityChange * Vector3(0, 100, 100), Vector4(0, 1, 0, 1), 0.0f);

		if (Window::GetKeyboard()->KeyPressed(NCL::KeyboardKeys::SPACE) && hookTimer > 50) {
			world->RemoveGameObject(Hook);
			Vector3 pos = (Ball->GetTransform()).GetPosition() + Vector3(((SphereVolume&) * (Ball->GetBoundingVolume())).GetRadius() + 0.1, 0, 0);

			Ray ray = Ray(pos, pos + gravityChange * Vector3(0, 1000, 1000));
			RayCollision closestCollision;
			if (world->Raycast(ray, closestCollision, true)) {
				Hook = AddSphereToWorld(closestCollision.collidedAt, Vector4(1, 1, 1, 1), 1, false, 0);
				hookConstraint = new PositionConstraint(Hook, Ball, (Hook->GetTransform().GetPosition() - Ball->GetTransform().GetPosition()).Length() - 20);
				world->AddConstraint(hookConstraint);
				hookTimer = 0;
			}
		}

		if (Window::GetKeyboard()->KeyDown(NCL::KeyboardKeys::SPACE) && Hook) {
			Debug::DrawLine(Ball->GetTransform().GetPosition(), Hook->GetTransform().GetPosition(), Vector4(0.4, 0.3, 0.3, 1), 0.0f);
		}

		if (!(Window::GetKeyboard()->KeyDown(NCL::KeyboardKeys::SPACE)) && Hook) {
			world->RemoveGameObject(Hook);
			world->RemoveConstraint(hookConstraint);
		}

		hookTimer++;
	}
}

/*

Every frame, this code will let you perform a raycast, to see if there's an object
underneath the cursor, and if so 'select it' into a pointer, so that it can be 
manipulated later. Pressing Q will let you toggle between this behaviour and instead
letting you move the camera around. 

*/
bool TutorialGame::SelectObject() {
	if (Window::GetKeyboard()->KeyPressed(KeyboardKeys::Q)) {
		inSelectionMode = !inSelectionMode;
		if (inSelectionMode) {
			Window::GetWindow()->ShowOSPointer(true);
			Window::GetWindow()->LockMouseToWindow(false);
		}
		else {
			Window::GetWindow()->ShowOSPointer(false);
			Window::GetWindow()->LockMouseToWindow(true);
		}
	}

	if (inSelectionMode) {

		if (Window::GetMouse()->ButtonDown(NCL::MouseButtons::LEFT)) {
			if (selectionObject) {	//set colour to deselected;
				selectionObject->GetRenderObject()->SetColour(Ocolour);
				if (zObject) {
					zObject->GetRenderObject()->SetColour(Vector4(1, 1, 1, 1));
					zObject = nullptr;
				}
				selectionObject = nullptr;
				lockedObject	= nullptr;
			}

			Ray ray = CollisionDetection::BuildRayFromMouse(*world->GetMainCamera());

			RayCollision closestCollision;
			RayCollision closestObject;
			if (world->Raycast(ray, closestCollision, true)) {
				selectionObject = (GameObject*)closestCollision.node;
				Ocolour = selectionObject->GetRenderObject()->GetColour();
				selectionObject->GetRenderObject()->SetColour(Vector4(0, 1, 0, 1));
				Vector3 pose = Vector3(0, 0, 0);
				return true;
			}
			else {
				return false;
			}
		}
	}

	if(selectionObject){
		GiveObjectInfo(selectionObject);
	}

	if (Window::GetKeyboard()->KeyPressed(NCL::KeyboardKeys::L)) {
		if (selectionObject) {
			if (lockedObject == selectionObject) {
				lockedObject = nullptr;
			}
			else {
				lockedObject = selectionObject;
			}
		}

	}

	return false;
}

/*
If an object has been clicked, it can be pushed with the right mouse button, by an amount
determined by the scroll wheel. In the first tutorial this won't do anything, as we haven't
added linear motion into our physics system. After the second tutorial, objects will move in a straight
line - after the third, they'll be able to twist under torque aswell.
*/
void TutorialGame::MoveSelectedObject() {
	//renderer->DrawString("Click Force:" + std::to_string(forceMagnitude),
	//Vector2(5, 50)); //Draw debug text at 10,20
	forceMagnitude += Window::GetMouse()->GetWheelMovement() * 100.0f;

	if (!selectionObject) {
		return;//we haven’t selected anything!
	}

	if (Window::GetMouse()->ButtonPressed(NCL::MouseButtons::LEFT)) {
		Ray ray = CollisionDetection::BuildRayFromMouse(*world->GetMainCamera());
		RayCollision closestCollision;
		if (world->Raycast(ray, closestCollision, true)) {
			if (closestCollision.node == selectionObject) {
				switch (selectionObject->GetType()) {
				case ObjectType::Push:
					if (selectionObject->GetRenderObject()->GetTransform()->GetPosition().y - 25 < Ball->GetRenderObject()->GetTransform()->GetPosition().y
						&& Ball->GetRenderObject()->GetTransform()->GetPosition().y < selectionObject->GetRenderObject()->GetTransform()->GetPosition().y + 25
						&& selectionObject->GetRenderObject()->GetTransform()->GetPosition().z - 10 < Ball->GetRenderObject()->GetTransform()->GetPosition().z
						&& Ball->GetRenderObject()->GetTransform()->GetPosition().z < selectionObject->GetRenderObject()->GetTransform()->GetPosition().z + 10)
					{
						Ball->GetPhysicsObject()->AddForce(Vector3(0, 2, 0) * forceMagnitude);
					}
						
					break;
				case ObjectType::Pull:
					if (selectionObject->GetRenderObject()->GetTransform()->GetPosition().y - 10 < Ball->GetRenderObject()->GetTransform()->GetPosition().y
						&& Ball->GetRenderObject()->GetTransform()->GetPosition().y < selectionObject->GetRenderObject()->GetTransform()->GetPosition().y + 10)
					{
						if(Ball->GetRenderObject()->GetTransform()->GetPosition().z > selectionObject->GetRenderObject()->GetTransform()->GetPosition().z)
							Ball->GetPhysicsObject()->AddForce(Vector3(0, 0, -2) * forceMagnitude);
						else
							Ball->GetPhysicsObject()->AddForce(Vector3(0, 0, 2) * forceMagnitude);
					}
					break;
				case ObjectType::Spin:
					selectionObject->GetPhysicsObject()->SetAngularVelocity(Vector3(10, 0, 0));
					break;
				case ObjectType::Gravity:
					changeGravity = !changeGravity;
					break;
				case ObjectType::SurfaceChanger:
					selectionObject->ChangeSurface();
					break;
				}	
			}
		}
	}


	//Push the selected object!
	if (Window::GetMouse()->ButtonPressed(NCL::MouseButtons::RIGHT)) {
		Ray ray = CollisionDetection::BuildRayFromMouse(* world->GetMainCamera());
		RayCollision closestCollision;
		if (world->Raycast(ray, closestCollision, true)) {
			if (closestCollision.node == selectionObject) {
				selectionObject->GetPhysicsObject()->AddForceAtPosition(ray.GetDirection() * forceMagnitude, closestCollision.collidedAt);
			}
		}
	}

	if (Window::GetKeyboard()->KeyPressed(NCL::KeyboardKeys::I)) {
		Ray ray = CollisionDetection::BuildRayFromMouse(*world->GetMainCamera());
		RayCollision closestCollision;
		if (world->Raycast(ray, closestCollision, true)) {
			if (closestCollision.node == selectionObject) {
				selectionObject->GetPhysicsObject()->AddForce(Vector3(1,0,0) * forceMagnitude);
			}
		}
	}
	if (Window::GetKeyboard()->KeyPressed(NCL::KeyboardKeys::J)) {
		Ray ray = CollisionDetection::BuildRayFromMouse(*world->GetMainCamera());
		RayCollision closestCollision;
		if (world->Raycast(ray, closestCollision, true)) {
			if (closestCollision.node == selectionObject) {
				selectionObject->GetPhysicsObject()->AddForce(Vector3(0, 1, 0) * forceMagnitude);
			}
		}
	}
	if (Window::GetKeyboard()->KeyPressed(NCL::KeyboardKeys::K)) {
		Ray ray = CollisionDetection::BuildRayFromMouse(*world->GetMainCamera());
		RayCollision closestCollision;
		if (world->Raycast(ray, closestCollision, true)) {
			if (closestCollision.node == selectionObject) {
				selectionObject->GetPhysicsObject()->AddForce(Vector3(0, 0, 1) * forceMagnitude);
			}
		}
	}
}

void TutorialGame::ReactingBlocks() {

}

void TutorialGame::GiveObjectInfo(GameObject* Object) {
	renderer->DrawString("Position: " + std::to_string(Object->GetRenderObject()->GetTransform()->GetPosition().x) + ", "
									  + std::to_string(Object->GetRenderObject()->GetTransform()->GetPosition().y) + ", "
									  + std::to_string(Object->GetRenderObject()->GetTransform()->GetPosition().z), Vector2(5, 5));
	renderer->DrawString("Orientation: " + std::to_string(Object->GetRenderObject()->GetTransform()->GetOrientation().x) + ", "
										 + std::to_string(Object->GetRenderObject()->GetTransform()->GetOrientation().y) + ", "
										 + std::to_string(Object->GetRenderObject()->GetTransform()->GetOrientation().z) + ", "
										 + std::to_string(Object->GetRenderObject()->GetTransform()->GetOrientation().w), Vector2(5, 10));
	renderer->DrawString("Elasticity: " + std::to_string(Object->GetPhysicsObject()->GetElasticity()), Vector2(5, 15));
	renderer->DrawString("Friction: " + std::to_string(Object->GetPhysicsObject()->GetFriction()), Vector2(5, 20));
	
	if(!GameB)
		renderer->DrawString("Object Type: " + Object->GetTypeName(), Vector2(5, 25));
}

void TutorialGame::RedBehaviourTree() {
	float behaviourTimer;
	float distanceToTarget;
	BehaviourAction* findKey = new BehaviourAction("Find Key",
		[&](float dt, BehaviourState state)->BehaviourState {
			if (state == Initialise) {
				target = TestPaths(Keys);
				state = Ongoing;
			}

			else if (state == Ongoing) {
				if(abs((target->GetRenderObject()->GetTransform()->GetPosition() - Red->GetRenderObject()->GetTransform()->GetPosition()).x) < 10 
					&& abs((target->GetRenderObject()->GetTransform()->GetPosition() - Red->GetRenderObject()->GetTransform()->GetPosition()).z) < 10){
					world->RemoveGameObject(target);
					Keys.erase(std::remove(Keys.begin(), Keys.end(), target), Keys.end());
					EnemyKeys++;
					return Success;
				}
				if (Pathfinding(Ball).size() < 10) {
					return Failure;
				}
			}
			return state; //will be ’ongoing ’ until success
		});
	BehaviourAction* findChest = new BehaviourAction("Find Chest",
		[&](float dt, BehaviourState state)->BehaviourState {
			if (state == Initialise) {
				target = TestPaths(Chests);
				state = Ongoing;
			}
			else if (state == Ongoing) {
				if (abs((target->GetRenderObject()->GetTransform()->GetPosition() - Red->GetRenderObject()->GetTransform()->GetPosition()).x) < 15
					&& abs((target->GetRenderObject()->GetTransform()->GetPosition() - Red->GetRenderObject()->GetTransform()->GetPosition()).z) < 15) {
					world->RemoveGameObject(target);
					Chests.erase(std::remove(Chests.begin(), Chests.end(), target), Chests.end());
					EnemyKeys--;
					EnemyPoints+=10;
					return Success;
				}
				if (Pathfinding(Ball).size() < 10) {
					return Failure;
				}
			}
			return state; //will be ’ongoing ’ until success
		});
	BehaviourAction* lookForPlayer = new BehaviourAction(
		"Look For Treasure",
		[&](float dt, BehaviourState state)->BehaviourState {
			if (state == Initialise) {
				target = Ball;
				return Ongoing;
			}
			else if (state == Ongoing) {
				if (Pathfinding(Ball).size() >= 10) {
					return Failure;
				}
				if (abs((target->GetRenderObject()->GetTransform()->GetPosition() - Red->GetRenderObject()->GetTransform()->GetPosition()).x) < 15
					&& abs((target->GetRenderObject()->GetTransform()->GetPosition() - Red->GetRenderObject()->GetTransform()->GetPosition()).z) < 15) {
					target->GetRenderObject()->GetTransform()->SetPosition(Vector3(20, 20, 20));
					Red->GetRenderObject()->GetTransform()->SetPosition(Vector3(380, 20, 380));
					PlayerPoints-=10;
					return Success;
				}
			}
			return state;
		});

	sequence->AddChild(findKey);
	sequence->AddChild(findChest);

	selection->AddChild(sequence);
	selection->AddChild(lookForPlayer);
}

vector<Vector3> TutorialGame::Pathfinding(GameObject* target) {
	vector <Vector3 > testNodes;

	NavigationGrid grid("TestGrid1.txt", Red->GetRenderObject()->GetTransform()->GetPosition());

	NavigationPath outPath;

	Vector3 startPos(Red->GetRenderObject()->GetTransform()->GetPosition());
	Vector3 endPos(target->GetRenderObject()->GetTransform()->GetPosition());

	bool found = grid.FindPath(startPos, endPos, outPath);

	Vector3 pos;
	while (outPath.PopWaypoint(pos)) {
		testNodes.push_back(pos);
	}
	
	return testNodes;
}

Vector3 TutorialGame::RedPathfinding(GameObject* target) {
	vector <Vector3 > testNodes = Pathfinding(target);

	if (testNodes.size() > 1)
		return testNodes[1];
	else if (testNodes.size() == 1)
		return testNodes[0];
	else
		return Vector3(0, 0, 0);
}

GameObject* TutorialGame::TestPaths(vector<GameObject*> targets) {
	float pathsize = FLT_MAX;
	GameObject* object = new GameObject();

	for (int i = 0; i < targets.size(); i++) {
		vector <Vector3 > testNodes = Pathfinding(targets[i]);

		if (testNodes.size() < pathsize) {
			pathsize = testNodes.size();
			object = targets[i];
		}
	}

	return object;
}

void TutorialGame::RedMove(GameObject* target) {
	if (abs((RedAim - Red->GetRenderObject()->GetTransform()->GetPosition()).x) < 1 && abs((RedAim - Red->GetRenderObject()->GetTransform()->GetPosition()).y) < 1 && abs((RedAim - Red->GetRenderObject()->GetTransform()->GetPosition()).z) < 1) {
		RedAim = RedPathfinding(target);
	}
	else
		Red->GetRenderObject()->GetTransform()->SetPosition(Red->GetRenderObject()->GetTransform()->GetPosition() + (RedAim - Red->GetRenderObject()->GetTransform()->GetPosition()).Normalised());
}
