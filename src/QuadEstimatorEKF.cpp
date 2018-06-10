#include "Common.h"
#include "QuadEstimatorEKF.h"
#include "Utility/SimpleConfig.h"
#include "Utility/StringUtils.h"
#include "Math/Quaternion.h"

using namespace SLR;

const int QuadEstimatorEKF::QUAD_EKF_NUM_STATES;

QuadEstimatorEKF::QuadEstimatorEKF(string config, string name)
  : BaseQuadEstimator(config),
  Q(QUAD_EKF_NUM_STATES, QUAD_EKF_NUM_STATES),
  R_GPS(6, 6),
  R_Mag(1, 1),
  ekfState(QUAD_EKF_NUM_STATES),
  ekfCov(QUAD_EKF_NUM_STATES, QUAD_EKF_NUM_STATES),
  trueError(QUAD_EKF_NUM_STATES)
{
	_name = name;
	Init();
}

QuadEstimatorEKF::~QuadEstimatorEKF()
{

}

void QuadEstimatorEKF::Init()
{
	ParamsHandle paramSys = SimpleConfig::GetInstance();

	paramSys->GetFloatVector(_config + ".InitState", ekfState);

	VectorXf initStdDevs(QUAD_EKF_NUM_STATES);
	paramSys->GetFloatVector(_config + ".InitStdDevs", initStdDevs);
	ekfCov.setIdentity();
	for (int i = 0; i < QUAD_EKF_NUM_STATES; i++)
	{
		ekfCov(i, i) = initStdDevs(i) * initStdDevs(i);
	}

	// complementary filter params
	attitudeTau = paramSys->Get(_config + ".AttitudeTau", .1f);
	dtIMU = paramSys->Get(_config + ".dtIMU", .002f);

	pitchEst = 0;
	rollEst = 0;
  
	// GPS measurement model covariance
	R_GPS.setZero();
	R_GPS(0, 0) = R_GPS(1, 1) = powf(paramSys->Get(_config + ".GPSPosXYStd", 0), 2);
	R_GPS(2, 2) = powf(paramSys->Get(_config + ".GPSPosZStd", 0), 2);
	R_GPS(3, 3) = R_GPS(4, 4) = powf(paramSys->Get(_config + ".GPSVelXYStd", 0), 2);
	R_GPS(5, 5) = powf(paramSys->Get(_config + ".GPSVelZStd", 0), 2);

	// magnetometer measurement model covariance
	R_Mag.setZero();
	R_Mag(0, 0) = powf(paramSys->Get(_config + ".MagYawStd", 0), 2);

	// load the transition model covariance
	Q.setZero();
	Q(0, 0) = Q(1, 1) = powf(paramSys->Get(_config + ".QPosXYStd", 0), 2);
	Q(2, 2) = powf(paramSys->Get(_config + ".QPosZStd", 0), 2);
	Q(3, 3) = Q(4, 4) = powf(paramSys->Get(_config + ".QVelXYStd", 0), 2);
	Q(5, 5) = powf(paramSys->Get(_config + ".QVelZStd", 0), 2);
	Q(6, 6) = powf(paramSys->Get(_config + ".QYawStd", 0), 2);
	Q *= dtIMU;

	rollErr = pitchErr = maxEuler = 0;
	posErrorMag = velErrorMag = 0;
}

void QuadEstimatorEKF::UpdateFromIMU(V3F accel, V3F gyro)
{
	// Improve a complementary filter-type attitude filter
	// 
	// Currently a small-angle approximation integration method is implemented
	// The integrated (predicted) value is then updated in a complementary filter style with attitude information from accelerometers
	// 
	// Implement a better integration method that uses the current attitude estimate (rollEst, pitchEst and ekfState(6))
	// to integrate the body rates into new Euler angles.
	//
	// HINTS:
	//  - there are several ways to go about this, including:
	//    1) create a rotation matrix based on your current Euler angles, integrate that, convert back to Euler angles
//    OR 
//    2) use the Quaternion<float> class, which has a handy FromEuler123_RPY function for creating a quaternion from Euler Roll/PitchYaw
//       (Quaternion<float> also has a IntegrateBodyRate function, though this uses quaternions, not Euler angles)

////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
// SMALL ANGLE GYRO INTEGRATION:
// (replace the code below)
// make sure you comment it out when you add your own code -- otherwise e.g. you might integrate yaw twice

// Create 3x3 matrix to convert body fixed angular rate measurements from gyro into the derivitive of Euler angles

float cosPhi = cos(rollEst);
float sinPhi = sin(rollEst);
float tanTheta = tan(pitchEst);
float secTheta = 1.0f / cos(pitchEst);

Mat3x3F R;
R(0, 0) = 1.0;
R(0, 1) = sinPhi * tanTheta;
R(0, 2) = cosPhi * tanTheta;
R(1, 0) = 0.0;
R(1, 1) = cosPhi;
R(1, 2) = -sinPhi;
R(2, 0) = 0.0;
R(2, 1) = sinPhi * secTheta;
R(2, 2) = cosPhi * secTheta;

V3F EulerRates = R * gyro;

float predictedPitch = pitchEst + dtIMU * EulerRates.y;
float predictedRoll = rollEst + dtIMU * EulerRates.x;
ekfState(6) += dtIMU * EulerRates.z;	// yaw

// since yaw is periodic force yaw state to stay in the range  -pi .. pi

// Note: the following line should do what the following lines do,.....but the fmodf function isn't working properly
//ekfState(6) = fmodf((ekfState(6) + M_PI), (2.0f * M_PI)) - M_PI;

float pi = float(M_PI);
if (abs(ekfState(6)) > pi)
{
	int cnt = 0;
	while ((cnt < 10) && (ekfState(6) < -pi))
	{
		ekfState(6) += 2.0f * pi;
		cnt++;
	}
	while ((cnt < 10) && (ekfState(6) > pi))
	{
		ekfState(6) -= 2.0f * pi;
		cnt++;
	}

	assert(abs(ekfState(6)) < M_PI);
}

/////////////////////////////// END STUDENT CODE ////////////////////////////

// CALCULATE UPDATE
accelRoll = atan2f(accel.y, accel.z);
accelPitch = atan2f(-accel.x, 9.81f);

// FUSE INTEGRATION AND UPDATE
rollEst = attitudeTau / (attitudeTau + dtIMU) * (predictedRoll)+dtIMU / (attitudeTau + dtIMU) * accelRoll;
pitchEst = attitudeTau / (attitudeTau + dtIMU) * (predictedPitch)+dtIMU / (attitudeTau + dtIMU) * accelPitch;

lastGyro = gyro;
}

void QuadEstimatorEKF::UpdateTrueError(V3F truePos, V3F trueVel, Quaternion<float> trueAtt)
{
	VectorXf trueState(QUAD_EKF_NUM_STATES);
	trueState(0) = truePos.x;
	trueState(1) = truePos.y;
	trueState(2) = truePos.z;
	trueState(3) = trueVel.x;
	trueState(4) = trueVel.y;
	trueState(5) = trueVel.z;
	trueState(6) = trueAtt.Yaw();

	trueError = ekfState - trueState;
	if (trueError(6) > F_PI) trueError(6) -= 2.f*F_PI;
	if (trueError(6) < -F_PI) trueError(6) += 2.f*F_PI;

	pitchErr = pitchEst - trueAtt.Pitch();
	rollErr = rollEst - trueAtt.Roll();
	maxEuler = MAX(fabs(pitchErr), MAX(fabs(rollErr), fabs(trueError(6))));

	posErrorMag = truePos.dist(V3F(ekfState(0), ekfState(1), ekfState(2)));
	velErrorMag = trueVel.dist(V3F(ekfState(3), ekfState(4), ekfState(5)));
}

VectorXf QuadEstimatorEKF::PredictState(VectorXf curState, float dt, V3F accel, V3F gyro)
{
	assert(curState.size() == QUAD_EKF_NUM_STATES);
	VectorXf predictedState = curState;
	// Predict the current state forward by time dt using current accelerations and body rates as input
	// INPUTS: 
	//   curState: starting state
	//   dt: time step to predict forward by [s]
	//   accel: acceleration of the vehicle, in body frame, *not including gravity* [m/s2]
	//   gyro: body rates of the vehicle, in body frame [rad/s]
	//   
	// OUTPUT:
	//   return the predicted state as a vector

	// HINTS 
	// - dt is the time duration for which you should predict. It will be very short (on the order of 1ms)
	//   so simplistic integration methods are fine here
	// - we've created an Attitude Quaternion for you from the current state. Use 
	//   attitude.Rotate_BtoI(<V3F>) to rotate a vector from body frame to inertial frame
	// - the yaw integral is already done in the IMU update. Be sure not to integrate it again here

	Quaternion<float> attitude = Quaternion<float>::FromEuler123_RPY(rollEst, pitchEst, curState(6));

	////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

	V3F accel_NED = attitude.Rotate_BtoI(accel) + V3F(0.0, 0.0, -9.81f);

	for (uint8_t k = 0; k < 3; k++)
	{
		// integrate velocity states to update position states
		predictedState[k] += dt * curState[k + 3];

		// integrate NED acceleration with gracity removed to update velocity states
		predictedState[k + 3] += dt * accel_NED[k];
	}

	/////////////////////////////// END STUDENT CODE ////////////////////////////

	return predictedState;
}

MatrixXf QuadEstimatorEKF::GetRbgPrime(float roll, float pitch, float yaw)
{
	// first, figure out the Rbg_prime
	MatrixXf RbgPrime(3, 3);
	RbgPrime.setZero();

	// Return the partial derivative of the Rbg rotation matrix with respect to yaw. We call this RbgPrime.
	// INPUTS: 
	//   roll, pitch, yaw: Euler angles at which to calculate RbgPrime
	//   
	// OUTPUT:
	//   return the 3x3 matrix representing the partial derivative at the given point

	// HINTS
	// - this is just a matter of putting the right sin() and cos() functions in the right place.
	//   make sure you write clear code and triple-check your math
	// - You can also do some numerical partial derivatives in a unit test scheme to check 
	//   that your calculations are reasonable

	////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

	float cphi = cos(rollEst);
	float sphi = sin(rollEst);
	float ctheta = cos(pitchEst);
	float stheta = sin(pitchEst);
	float cpsi = cos(ekfState[6]);
	float spsi = sin(ekfState[6]);

	RbgPrime(0, 0) = -ctheta * spsi;
	RbgPrime(1, 0) = ctheta * cpsi;

	RbgPrime(0, 0) = -sphi * stheta * spsi - cphi * cpsi;
	RbgPrime(1, 0) = sphi * stheta * cpsi - cphi * spsi;

	RbgPrime(0, 0) = -cphi * stheta * spsi + sphi * cpsi;
	RbgPrime(1, 0) = cphi * stheta * cpsi + sphi * spsi;

	/////////////////////////////// END STUDENT CODE ////////////////////////////

	return RbgPrime;
}

void QuadEstimatorEKF::Predict(float dt, V3F accel, V3F gyro)
{
	// predict the state forward
	VectorXf newState = PredictState(ekfState, dt, accel, gyro);

	// Predict the current covariance forward by dt using the current accelerations and body rates as input.
	// INPUTS: 
	//   dt: time step to predict forward by [s]
	//   accel: acceleration of the vehicle, in body frame, *not including gravity* [m/s2]
	//   gyro: body rates of the vehicle, in body frame [rad/s]
	//   state (member variable): current state (state at the beginning of this prediction)
	//   
	// OUTPUT:
	//   update the member variable cov to the predicted covariance

	// HINTS
	// - update the covariance matrix cov according to the EKF equation.
	// 
	// - you may find the current estimated attitude in variables rollEst, pitchEst, state(6).
	//
	// - use the class MatrixXf for matrices. To create a 3x5 matrix A, use MatrixXf A(3,5).
	//
	// - the transition model covariance, Q, is loaded up from a parameter file in member variable Q
	// 
	// - This is unfortunately a messy step. Try to split this up into clear, manageable steps:
	//   1) Calculate the necessary helper matrices, building up the transition jacobian
	//   2) Once all the matrices are there, write the equation to update cov.
	//
	// - if you want to transpose a matrix in-place, use A.transposeInPlace(), not A = A.transpose()
	// 

	// we'll want the partial derivative of the Rbg matrix
	MatrixXf RbgPrime = GetRbgPrime(rollEst, pitchEst, ekfState(6));

	// we've created an empty Jacobian for you, currently simply set to identity
	MatrixXf gPrime(QUAD_EKF_NUM_STATES, QUAD_EKF_NUM_STATES);
	gPrime.setIdentity();

	////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

	MatrixXf tmp(3, 1);
	tmp(0, 0) = accel.x;
	tmp(1, 0) = accel.y;
	tmp(2, 0) = accel.z;
	tmp = RbgPrime * tmp;

	for (uint8_t k = 0; k < 7; k++)
	{
		gPrime(k, k) = 1.0;

		if (k < 3)
		{
			gPrime(k, k + 3) = dt;
			gPrime(k + 3, 6) = tmp(k) * dt;
		}
	}

	// update covariance from  prediction
	ekfCov = gPrime * ekfCov * gPrime.transpose() + Q;

	/////////////////////////////// END STUDENT CODE ////////////////////////////

	ekfState = newState;
}

void QuadEstimatorEKF::UpdateFromGPS(V3F pos, V3F vel)
{
	VectorXf z(6), zFromX(6);
	z(0) = pos.x;
	z(1) = pos.y;
	z(2) = pos.z;
	z(3) = vel.x;
	z(4) = vel.y;
	z(5) = vel.z;

	MatrixXf hPrime(6, QUAD_EKF_NUM_STATES);
	hPrime.setZero();

	// GPS UPDATE
	// Hints: 
	//  - The GPS measurement covariance is available in member variable R_GPS
	//  - this is a very simple update
	////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

	/////////////////////////////// END STUDENT CODE ////////////////////////////

	Update(z, hPrime, R_GPS, zFromX);
}

void QuadEstimatorEKF::UpdateFromMag(float magYaw)
{
	VectorXf z(1), zFromX(1);
	z(0) = magYaw;

	MatrixXf hPrime(1, QUAD_EKF_NUM_STATES);
	hPrime.setZero();

	float pi = float(M_PI);

	// MAGNETOMETER UPDATE
	// Hints: 
	//  - Your current estimated yaw can be found in the state vector: ekfState(6)
	//  - Make sure to normalize the difference between your measured and estimated yaw
	//    (you don't want to update your yaw the long way around the circle)
	//  - The magnetomer measurement covariance is available in member variable R_Mag
	////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

	hPrime(6) = 1.0;
	
	
	// setting zFromX(0) to zero effectively makes z the error and not the measurement
	zFromX(0) = ekfState(6);

	// The following lines should be able to be replace with a single lines using fmodf, but fmodf if broken in the environment
	
	if (abs(z(0) - zFromX(0)) > pi)
	{
		int cnt = 0;
		float error = z(0) - zFromX(0);
		while ((cnt < 10) && (error < -pi))
		{
			error += 2.0f * pi;
			cnt++;
		}
		while ((cnt < 10) && (error > pi))
		{
			error -= 2.0f * pi;
			cnt++;
		}

		z(0) = zFromX(0) + error;

		assert(abs(zFromX(0) - z(0)) < M_PI);
	}

	/////////////////////////////// END STUDENT CODE ////////////////////////////

	Update(z, hPrime, R_Mag, zFromX);

	// since yaw is periodic force yaw state to stay in the range  -pi .. pi
	// The following lines should be able to be replace with a single lines using fmodf, but fmodf if broken in the environment
	if (abs(ekfState(6)) > pi)
	{
		int cnt = 0;
		while ((cnt < 10) && (ekfState(6) < -pi))
		{
			ekfState(6) += 2.0f * pi;
			cnt++;
		}
		while ((cnt < 10) && (ekfState(6) > pi))
		{
			ekfState(6) -= 2.0f * pi;
			cnt++;
		}

		assert(abs(ekfState(6)) < M_PI);
	}
}

// Execute an EKF update step
// z: measurement
// H: Jacobian of observation function evaluated at the current estimated state
// R: observation error model covariance 
// zFromX: measurement prediction based on current state
void QuadEstimatorEKF::Update(VectorXf& z, MatrixXf& H, MatrixXf& R, VectorXf& zFromX)
{
	assert(z.size() == H.rows());
	assert(QUAD_EKF_NUM_STATES == H.cols());
	assert(z.size() == R.rows());
	assert(z.size() == R.cols());
	assert(z.size() == zFromX.size());

	MatrixXf toInvert(z.size(), z.size());
	toInvert = H*ekfCov*H.transpose() + R;
	MatrixXf K = ekfCov * H.transpose() * toInvert.inverse();

	ekfState = ekfState + K*(z - zFromX);

	MatrixXf eye(QUAD_EKF_NUM_STATES, QUAD_EKF_NUM_STATES);
	eye.setIdentity();

	ekfCov = (eye - K*H)*ekfCov;
}

// Calculate the condition number of the EKF ovariance matrix (useful for numerical diagnostics)
// The condition number provides a measure of how similar the magnitudes of the error metric beliefs 
// about the different states are. If the magnitudes are very far apart, numerical issues will start to come up.
float QuadEstimatorEKF::CovConditionNumber() const
{
	MatrixXf m(7, 7);
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			m(i, j) = ekfCov(i, j);
		}
	}

	Eigen::JacobiSVD<MatrixXf> svd(m);
	float cond = svd.singularValues()(0)
	/ svd.singularValues()(svd.singularValues().size() - 1);
	return cond;
}

// Access functions for graphing variables
bool QuadEstimatorEKF::GetData(const string& name, float& ret) const
{
	if (name.find_first_of(".") == string::npos) return false;
	string leftPart = LeftOf(name, '.');
	string rightPart = RightOf(name, '.');

	if (ToUpper(leftPart) == ToUpper(_name))
	{
#define GETTER_HELPER(A,B) if (SLR::ToUpper(rightPart) == SLR::ToUpper(A)){ ret=(B); return true; }
		GETTER_HELPER("Est.roll", rollEst);
		GETTER_HELPER("Est.pitch", pitchEst);

		GETTER_HELPER("Est.x", ekfState(0));
		GETTER_HELPER("Est.y", ekfState(1));
		GETTER_HELPER("Est.z", ekfState(2));
		GETTER_HELPER("Est.vx", ekfState(3));
		GETTER_HELPER("Est.vy", ekfState(4));
		GETTER_HELPER("Est.vz", ekfState(5));
		GETTER_HELPER("Est.yaw", ekfState(6));

		GETTER_HELPER("Est.S.x", sqrtf(ekfCov(0, 0)));
		GETTER_HELPER("Est.S.y", sqrtf(ekfCov(1, 1)));
		GETTER_HELPER("Est.S.z", sqrtf(ekfCov(2, 2)));
		GETTER_HELPER("Est.S.vx", sqrtf(ekfCov(3, 3)));
		GETTER_HELPER("Est.S.vy", sqrtf(ekfCov(4, 4)));
		GETTER_HELPER("Est.S.vz", sqrtf(ekfCov(5, 5)));
		GETTER_HELPER("Est.S.yaw", sqrtf(ekfCov(6, 6)));

		// diagnostic variables
		GETTER_HELPER("Est.D.AccelPitch", accelPitch);
		GETTER_HELPER("Est.D.AccelRoll", accelRoll);

		GETTER_HELPER("Est.D.ax_g", accelG[0]);
		GETTER_HELPER("Est.D.ay_g", accelG[1]);
		GETTER_HELPER("Est.D.az_g", accelG[2]);

		GETTER_HELPER("Est.E.x", trueError(0));
		GETTER_HELPER("Est.E.y", trueError(1));
		GETTER_HELPER("Est.E.z", trueError(2));
		GETTER_HELPER("Est.E.vx", trueError(3));
		GETTER_HELPER("Est.E.vy", trueError(4));
		GETTER_HELPER("Est.E.vz", trueError(5));
		GETTER_HELPER("Est.E.yaw", trueError(6));
		GETTER_HELPER("Est.E.pitch", pitchErr);
		GETTER_HELPER("Est.E.roll", rollErr);
		GETTER_HELPER("Est.E.MaxEuler", maxEuler);

		GETTER_HELPER("Est.E.pos", posErrorMag);
		GETTER_HELPER("Est.E.vel", velErrorMag);

		GETTER_HELPER("Est.D.covCond", CovConditionNumber());
	#undef GETTER_HELPER
	}
	return false;
};

vector<string> QuadEstimatorEKF::GetFields() const
{
  vector<string> ret = BaseQuadEstimator::GetFields();
  ret.push_back(_name + ".Est.roll");
  ret.push_back(_name + ".Est.pitch");

  ret.push_back(_name + ".Est.x");
  ret.push_back(_name + ".Est.y");
  ret.push_back(_name + ".Est.z");
  ret.push_back(_name + ".Est.vx");
  ret.push_back(_name + ".Est.vy");
  ret.push_back(_name + ".Est.vz");
  ret.push_back(_name + ".Est.yaw");

  ret.push_back(_name + ".Est.S.x");
  ret.push_back(_name + ".Est.S.y");
  ret.push_back(_name + ".Est.S.z");
  ret.push_back(_name + ".Est.S.vx");
  ret.push_back(_name + ".Est.S.vy");
  ret.push_back(_name + ".Est.S.vz");
  ret.push_back(_name + ".Est.S.yaw");

  ret.push_back(_name + ".Est.E.x");
  ret.push_back(_name + ".Est.E.y");
  ret.push_back(_name + ".Est.E.z");
  ret.push_back(_name + ".Est.E.vx");
  ret.push_back(_name + ".Est.E.vy");
  ret.push_back(_name + ".Est.E.vz");
  ret.push_back(_name + ".Est.E.yaw");
  ret.push_back(_name + ".Est.E.pitch");
  ret.push_back(_name + ".Est.E.roll");

  ret.push_back(_name + ".Est.E.pos");
  ret.push_back(_name + ".Est.E.vel");

  ret.push_back(_name + ".Est.E.maxEuler");

  ret.push_back(_name + ".Est.D.covCond");

  // diagnostic variables
  ret.push_back(_name + ".Est.D.AccelPitch");
  ret.push_back(_name + ".Est.D.AccelRoll");
  ret.push_back(_name + ".Est.D.ax_g");
  ret.push_back(_name + ".Est.D.ay_g");
  ret.push_back(_name + ".Est.D.az_g");
  return ret;
};
