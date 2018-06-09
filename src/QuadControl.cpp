#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

void QuadControl::Init()
{
	BaseController::Init();

	// variables needed for integral control
	integratedAltitudeError = 0;

#ifndef __PX4_NUTTX
	// Load params from simulator parameter system
	ParamsHandle config = SimpleConfig::GetInstance();

	// Load parameters (default to 0)
	kpPosXY = config->Get(_config + ".kpPosXY", 0);
	kpPosZ = config->Get(_config + ".kpPosZ", 0);
	KiPosZ = config->Get(_config + ".KiPosZ", 0);

	kpVelXY = config->Get(_config + ".kpVelXY", 0);
	kpVelZ = config->Get(_config + ".kpVelZ", 0);

	kpBank = config->Get(_config + ".kpBank", 0);
	kpYaw = config->Get(_config + ".kpYaw", 0);

	kpPQR = config->Get(_config + ".kpPQR", V3F());

	maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
	maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
	maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
	maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

	maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

	minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
	maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
	// load params from PX4 parameter system
	//TODO
	param_get(param_find("MC_PITCH_P"), &Kp_bank);
	param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
	// Convert a desired 3-axis moment and collective thrust command to 
	//   individual motor thrust commands
	// INPUTS: 
	//   desCollectiveThrust: desired collective thrust [N]
	//   desMoment: desired rotation moment about each axis [N m]
	// OUTPUT:
	//   set class member variable cmd (class variable for graphing) where
	//   cmd.desiredThrustsN[0..3]: motor commands, in [N]

	// HINTS: 
	// - you can access parts of desMoment via e.g. desMoment.x
	// You'll need the arm length parameter L, and the drag/thrust ratio kappa

 	////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

	float L_eff = L / sqrt(2.0f);  // effive length or moment arm to blade in both x & y
	float Fc_blade = collThrustCmd / 4.0f; // average for per blade to generate collective thrust
	float dF_Mx_blade = momentCmd.x / L_eff / 4.0f; // delta force per blade to create x moment
	float dF_My_blade = momentCmd.y / L_eff / 4.0f; // delta force per blade to create y moment
	float dF_Mz_blade = momentCmd.z / kappa / 4.0f; // delta force per blade to create z moment

	cmd.desiredThrustsN[0] = Fc_blade + dF_Mx_blade + dF_My_blade - dF_Mz_blade; // front left
	cmd.desiredThrustsN[1] = Fc_blade - dF_Mx_blade + dF_My_blade + dF_Mz_blade; // front right
	cmd.desiredThrustsN[2] = Fc_blade + dF_Mx_blade - dF_My_blade + dF_Mz_blade; // rear left
	cmd.desiredThrustsN[3] = Fc_blade - dF_Mx_blade - dF_My_blade - dF_Mz_blade; // rear right



	/////////////////////////////// END STUDENT CODE ////////////////////////////
	
	return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
	// Calculate a desired 3-axis moment given a desired and current body rate
	// INPUTS: 
	//   pqrCmd: desired body rates [rad/s]
	//   pqr: current or estimated body rates [rad/s]
	// OUTPUT:
	//   return a V3F containing the desired moments for each of the 3 axes

	// HINTS: 
	//  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
	//  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
	//  - you'll also need the gain parameter kpPQR (it's a V3F)

	// use p gain to calculate pqrcmd_dot or the desired rotation acceleration rates
	
	V3F momentCmd;

	////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

	V3F pqrCmd_dot = kpPQR * (pqrCmd - pqr);

	momentCmd.x = Ixx * pqrCmd_dot.x + (Izz - Iyy) * pqr.z *  pqr.y;
	momentCmd.y = Iyy * pqrCmd_dot.y + (Ixx - Izz) * pqr.x * pqr.z;
	momentCmd.z = Izz * pqrCmd_dot.z + (Iyy - Ixx) *  pqr.y * pqr.x;

	/////////////////////////////// END STUDENT CODE ////////////////////////////

	return momentCmd;
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
	// Calculate a desired pitch and roll angle rates based on a desired global
	//   lateral acceleration, the current attitude of the quad, and desired
	//   collective thrust command
	// INPUTS: 
	//   accelCmd: desired acceleration in global XY coordinates [m/s2]
	//   attitude: current or estimated attitude of the vehicle
	//   collThrustCmd: desired collective thrust of the quad [N]
	// OUTPUT:
	//   return a V3F containing the desired pitch and roll rates. The Z
	//     element of the V3F should be left at its default value (0)

	// HINTS: 
	//  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
	//  - you'll need the roll/pitch gain kpBank
	//  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first

	V3F pqrCmd;
	Mat3x3F R = attitude.RotationMatrix_IwrtB();

	////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

	// b actual x,y 
	float bx_a = R(0, 2); // R13
	float by_a = R(1, 2); // R23

	// b command x,y
	float bx_c = -CONSTRAIN(accelCmd.x * mass / max(minMotorThrust, collThrustCmd), -maxTiltAngle, maxTiltAngle);
	float by_c = -CONSTRAIN(accelCmd.y * mass / max(minMotorThrust, collThrustCmd), -maxTiltAngle, maxTiltAngle);

	// b command dot x,y
	float bx_c_dot = kpBank * (bx_c - bx_a);
	float by_c_dot = kpBank * (by_c - by_a);

	float R11 = R(0, 0);
	float R12 = R(0, 1);
	float R21 = R(1, 0);
	float R22 = R(1, 1);
	float R33 = R(2, 2);

	pqrCmd.x = (R21 * bx_c_dot - R11 * by_c_dot) / R33;
	pqrCmd.y = (R22 * bx_c_dot - R12 * by_c_dot) / R33;
	pqrCmd.z = 0.0f;

	/////////////////////////////// END STUDENT CODE ////////////////////////////

	return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
	// Calculate desired quad thrust based on altitude setpoint, actual altitude,
	//   vertical velocity setpoint, actual vertical velocity, and a vertical 
	//   acceleration feed-forward command
	// INPUTS: 
	//   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
	//   posZ, velZ: current vertical position and velocity in NED [m]
	//   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
	//   dt: the time step of the measurements [seconds]
	// OUTPUT:
	//   return a collective thrust command in [N]

	// HINTS: 
	//  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
	//  - you'll need the gain parameters kpPosZ and kpVelZ
	//  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
	//  - make sure to return a force, not an acceleration
	//  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

	Mat3x3F R = attitude.RotationMatrix_IwrtB();
	float thrust = 0;

	////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

	// position errors
	float dp = posZCmd - posZ;

	// integrator sum
	integratedAltitudeError += dp * dt;

	// add p control to velcoity command before limiting velocity....
	//   this means the effeictive p gain is the production of the pos and vel gain
	velZCmd += kpPosZ * dp; 
	float dv = CONSTRAIN(velZCmd, -maxAscentRate, maxDescentRate) - velZ;

	// desired "acceleration" command
	float u = accelZCmd + kpVelZ * dv + KiPosZ * integratedAltitudeError;

	// convert to total force command and divide by bz so the vertical component matches desired accel
	thrust = mass * (9.81f - u) / R(2, 2);

	/////////////////////////////// END STUDENT CODE ////////////////////////////

	return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmd)
{
	// Calculate a desired horizontal acceleration based on 
	//  desired lateral position/velocity/acceleration and current pose
	// INPUTS: 
	//   posCmd: desired position, in NED [m]
	//   velCmd: desired velocity, in NED [m/s]
	//   pos: current position, NED [m]
	//   vel: current velocity, NED [m/s]
	//   accelCmd: desired acceleration, NED [m/s2]
	// OUTPUT:
	//   return a V3F with desired horizontal accelerations. 
	//     the Z component should be 0
	// HINTS: 
	//  - use fmodf(foo,b) to constrain float foo to range [0,b]
	//  - use the gain parameters kpPosXY and kpVelXY
	//  - make sure you cap the horizontal velocity and acceleration
	//    to maxSpeedXY and maxAccelXY

	// make sure we don't have any incoming z-component
	accelCmd.z = 0;
	velCmd.z = 0;
	posCmd.z = pos.z;

	////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

	// add p control to velcoity command before limiting velocity....
	//   this means the effeictive p gain is the production of the pos and vel gain
	velCmd += kpPosXY * (posCmd - pos);

	// limit lateral velocity
	if (velCmd.mag() > maxSpeedXY)
	{
		velCmd *= maxSpeedXY / velCmd.mag();
	}

	// add D control 
	accelCmd += kpVelXY * (velCmd - vel);

	// limit lateral acceleration
	if (accelCmd.mag() > maxAccelXY)
	{
		accelCmd *= maxAccelXY / accelCmd.mag();
	}

	/////////////////////////////// END STUDENT CODE ////////////////////////////

	return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
	// Calculate a desired yaw rate to control yaw to yawCmd
	// INPUTS: 
	//   yawCmd: commanded yaw [rad]
	//   yaw: current yaw [rad]
	// OUTPUT:
	//   return a desired yaw rate [rad/s]
	// HINTS: 
	//  - use fmodf(foo,b) to constrain float foo to range [0,b] 
	//  - use the yaw control gain parameter kpYaw

	float yawRateCmd=0;

	////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

	yawRateCmd = kpYaw * (fmodf(((yawCmd - yaw) + M_PI), (2.0f * M_PI)) - M_PI);

	/////////////////////////////// END STUDENT CODE ////////////////////////////
  
	return yawRateCmd;
}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
	curTrajPoint = GetNextTrajectoryPoint(simTime);

	float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

	// reserve some thrust margin for angle control
	float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
	collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  
	V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
	V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
	desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

	V3F desMoment = BodyRateControl(desOmega, estOmega);

	return GenerateMotorCommands(collThrustCmd, desMoment);
}
