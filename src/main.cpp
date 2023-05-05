#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "rotation.h"
#include "earth.h"

using namespace std;
using namespace Eigen;

// record PVA information
struct PVA {
	Vector3d att;     // attitude in NED FRD
	Vector3d vel;     // velocity in NED 
	Vector3d pos;     // position in BLH
	Matrix3d Cbn;     // DCM downb_upn
	Quaterniond qbn;  // quaternion downb_upn
};


/**
	 * @brief INS Mechanization, update velocity using imudata
	 * @param [in]     pvapre	PVA information at the previous time (can be denoted as m - 1)
     *                       
	 * @param [in,out] pvacur   PVA information at the current time (can be denoted as m)
	 *                        
	 * @param [in]     imupre, imucur   IMU data at the previous time and at the current time
	 * 
	 * @param [in]     ts       update period
	 * */
void velUpdate(const PVA& pvapre, PVA& pvacur, const Vector<double, 7>& imupre, 
	           const Vector<double, 7>& imucur, const double& ts);


/**
	 * @brief INS Mechanization, update position using imudata
	 * @param [in]     pvapre	PVA information at the previous time (can be denoted as m - 1)
	 *
	 * @param [in,out] pvacur   PVA information at the current time (can be denoted as m)
	 *
	 * @param [in]     imupre, imucur   IMU data at the previous time and at the current time
	 *
	 * @param [in]     ts       update period
	 * */
void posUpdate(const PVA& pvapre, PVA& pvacur, const Vector<double, 7>& imupre,
	const Vector<double, 7>& imucur, const double& ts);


/**
	 * @brief INS Mechanization, update attitude using imudata
	 * @param [in]     pvapre	PVA information at the previous time (can be denoted as m - 1)
	 *
	 * @param [in,out] pvacur   PVA information at the current time (can be denoted as m)
	 *
	 * @param [in]     imupre, imucur   IMU data at the previous time and at the current time
	 *
	 * @param [in]     ts       update period
	 * */
void attUpdate(const PVA& pvapre, PVA& pvacur, const Vector<double, 7>& imupre,
	const Vector<double, 7>& imucur, const double& ts);


int main() {

	// 1. Load IMU data file, in .bin format 
	string imuname = "../data/imu.bin";   // .. represent the parent path
	ifstream imufile(imuname, ios_base::in | ios::binary);
	if (!imufile.is_open()) {
		cout << "The bin file of IMUData is Failed to load!\n";
		return -1;
	}	

	// 2. find the starting time because the imudata starts earlier
	double start_time = 91620.0;    // starting time (for this homework, is 91620.0 )
	Vector<double, 7> imudata;
	Vector<double, 7> imupre;
	while (imufile.read((char*)&imudata, sizeof(Vector<double, 7>))) {
		if (abs(imudata[0] - start_time) < 0.5e-3) {  // find the starting point successfully， do not using '==' for double format data
			imupre = imudata;
			break;
		}
	}
	if (imufile.eof()) {
		cout << "Failed to find the starting point! Please Check!\n";
		return -2;
	}


	// 3. set initial PVA
	double ts = 0.005;	     // sampling period (200Hz)
	PVA avp0;                // initial pva information
	Vector3d att0(D2R(0.0107951084511778), D2R(-2.14251290749072), D2R(-75.7498049314083)); // initial angle in order ZYX (NED:FRD)
	Vector3d vel0(0, 0, 0);
	Vector3d pos0(D2R(23.1373950708), D2R(113.3713651222), 2.175);  // initial pos in BLH
	avp0.att = att0;  avp0.vel = vel0; avp0.pos = pos0;  
	avp0.Cbn = Rotation::euler2matrix(att0);
	avp0.qbn = Rotation::euler2quaternion(att0);

	// 4. INS calculation using single simple + pre period method
	uint8_t n = 1;						          // the number of simple you chosen
	PVA avp;
	string res = "./result.txt";				  // write result to a txt file and then can using MATLAB for futher research easily
	ofstream resfile(res, ios_base::out);
	if (!resfile.is_open()) {
		cout << "Result File Failed to open!\n";
		return -3;
	}
	for (unsigned int i = 0; i < 1000; i += n) {  // 700,000+ epoches for INS calculation, now only choose first 1000 epoches
		cout << "Processing Line: " << i + 1 << endl;

		if (imufile.read((char*)&imudata, sizeof(Vector<double, 7>))) { // imucur
			// LOOK OUT: The Order cannot be changed
			velUpdate(avp0, avp, imupre, imudata, ts); // velocity update
			posUpdate(avp0, avp, imupre, imudata, ts); // position update
			attUpdate(avp0, avp, imupre, imudata, ts); // attitude update
			avp0 = avp;
			imupre = imudata;

			resfile << setprecision(15) << start_time + (i+1) * ts << "  ";
			resfile << setprecision(15) << R2D(avp.pos[0]) << "  " << R2D(avp.pos[1]) << "  " << avp.pos[2] << "  ";
			resfile << setprecision(15) << (avp.vel[0]) << "  " << (avp.vel[1]) << "  " << avp.vel[2] << "  ";
			resfile << setprecision(15) << R2D(avp.att[0]) << "  " << R2D(avp.att[1]) << "  " << R2D(avp.att[2]) << endl;
		}
		else
			break;
	}

	resfile.close();
	return 0;
}


void velUpdate(const PVA& pvapre, PVA& pvacur, const Vector<double, 7>& imupre,
	const Vector<double, 7>& imucur, const double& ts) {
	Vector3d d_vfb, d_vfn, d_vgn, gl, midvel, midpos;
	Vector3d temp1, temp2, temp3;
	Matrix3d cnn, I33 = Matrix3d::Identity();
	Quaterniond qne, qee, qnn, qbb, q1, q2;

	// 1. calculate geographic parameters, Meridian and Mao unitary radii,
    //    earth rotational angular velocity projected to n-frame,
    //    rotational angular velocity of n-frame to e-frame projected to n-frame, and gravity
	Vector2d rmrn = Earth::meridianPrimeVerticalRadius(pvapre.pos(0));
	Vector3d wie_n, wen_n;                                                            // NOTE: [Gongmin Yan's book ---- ENU]
	wie_n << WGS84_WIE * cos(pvapre.pos[0]), 0, -WGS84_WIE * sin(pvapre.pos[0]);         // Gongmin Yan (4.2.31a)
	wen_n << pvapre.vel[1] / (rmrn[1] + pvapre.pos[2]), -pvapre.vel[0] / (rmrn[0] + pvapre.pos[2]),
		-pvapre.vel[1] * tan(pvapre.pos[0]) / (rmrn[1] + pvapre.pos[2]);                 // Gongmin Yan (4.2.31b)
	double gravity = Earth::gravity(pvapre.pos);

	// 2. rotational and sculling motion
//	temp1 = imucur.dtheta.cross(imucur.dvel) / 2;
//	temp2 = imupre.dtheta.cross(imucur.dvel) / 12;
//	temp3 = imupre.dvel.cross(imucur.dtheta) / 12;
	temp1 = imucur.segment<3>(1).cross(imucur.segment<3>(4)) / 2;    // Gongmin Yan (4.1.36)
	temp2 = imupre.segment<3>(1).cross(imucur.segment<3>(4)) / 12;   // Gongmin Yan (4.1.44)
	temp3 = imupre.segment<3>(4).cross(imucur.segment<3>(1)) / 12;   // NOTE: [Gongmin Yan's book ---- two-samples method]
	// velocity increment due to the specific force in b-frame
//	d_vfb = imucur.dvel + temp1 + temp2 + temp3;
	d_vfb = imucur.segment<3>(4) + temp1 + temp2 + temp3;            // Gongmin Yan (4.1.32)

	// 3. velocity increment dut to the specfic force projected to the n-frame
	temp1 = (wie_n + wen_n) * ts / 2;
	cnn = I33 - Rotation::skewSymmetric(temp1);  
	d_vfn = cnn * pvapre.Cbn * d_vfb;								 // Gongmin Yan (4.1.50)

	// 4. velocity increment due to the gravity and Coriolis force
	gl << 0, 0, gravity;
	d_vgn = (gl - (2 * wie_n + wen_n).cross(pvapre.vel)) * ts;       // Gongmin Yan (4.1.27) 

	// 5. velocity at m-1/2
	midvel = pvapre.vel + (d_vfn + d_vgn) / 2;   // maybe this procedure is useless, I'll try to neglect it later
	// NOTE：pvapre.vel is in frame n(m-1), but d_vfn and d_vgn are in frame n(m); the difference between frame n(m-1) and n(m) are neglected!

	// position extrapolation to m-1/2
	qnn = Rotation::rotvec2quaternion(temp1);    // q_downn(m-1/2)_upn(m-1)
	temp2 << 0, 0, -WGS84_WIE * ts / 2;
	qee = Rotation::rotvec2quaternion(temp2);    // q_downe(m-1)_upe(m-1/2)
	qne = Earth::qne(pvapre.pos);                // q_downn(m-1)_upe(m-1)
	qne = qee * qne * qnn;                       // q_downn(m-1/2)_upe(m-1/2)
	midpos[2] = pvapre.pos[2] - midvel[2] * ts / 2;   // using '-' not '+' because of using FR'D'
	midpos = Earth::blh(qne, midpos[2]);

	// 6. recompute rmrn, wie_e, wen_n at m-1/2
	rmrn = Earth::meridianPrimeVerticalRadius(midpos[0]);
	wie_n << WGS84_WIE * cos(midpos[0]), 0, -WGS84_WIE * sin(midpos[0]);
	wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
		-midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);
	// recompute d_vfn
	temp3 = (wie_n + wen_n) * ts / 2;
	cnn = I33 - Rotation::skewSymmetric(temp3);
	d_vfn = cnn * pvapre.Cbn * d_vfb;
	// recompute d_vgn
	gl << 0, 0, Earth::gravity(midpos);
	d_vgn = (gl - (2 * wie_n + wen_n).cross(midvel)) * ts;

	// 7. velocity update finish
	pvacur.vel = pvapre.vel + d_vfn + d_vgn;
}


void posUpdate(const PVA& pvapre, PVA& pvacur, const Vector<double, 7>& imupre,
	const Vector<double, 7>& imucur, const double& ts) {

	Vector3d temp1, temp2, midvel, midpos;
	Quaterniond qne, qee, qnn;

	// 1. recompute velocity and position at m-1/2
	midvel = (pvacur.vel + pvapre.vel) / 2;
	midpos = pvapre.pos + Earth::DRi(pvapre.pos) * midvel * ts / 2;  // Gongmin Yan (4.1.57)

	// 2. recompute rmrn, wie_n, wen_n at m-1/2
	Vector2d rmrn;
	Vector3d wie_n, wen_n;
	rmrn = Earth::meridianPrimeVerticalRadius(midpos[0]);
	wie_n << WGS84_WIE * cos(midpos[0]), 0, -WGS84_WIE * sin(midpos[0]);
	wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
		-midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

	// 3. recompute n-frame rotation vector (n(m) with respect to n(m-1)-frame)
	temp1 = (wie_n + wen_n) * ts;
	qnn = Rotation::rotvec2quaternion(temp1);
	// e-frame rotation vector (e(m-1) with respect to e(m)-frame)
	temp2 << 0, 0, -WGS84_WIE * ts;             
	qee = Rotation::rotvec2quaternion(temp2);   // q_downe(m-1)_upe(m)

	// 4. position update finish
	qne = Earth::qne(pvapre.pos);
	qne = qee * qne * qnn;
	pvacur.pos[2] = pvapre.pos[2] - midvel[2] * ts;
	pvacur.pos = Earth::blh(qne, pvacur.pos[2]);
}


void attUpdate(const PVA& pvapre, PVA& pvacur, const Vector<double, 7>& imupre,
	const Vector<double, 7>& imucur, const double& ts) {
	Quaterniond qne_pre, qne_cur, qne_mid, qnn, qbb;
	Vector3d temp1, midpos, midvel;

	// 1. recompute velocity and position at m-1/2
	midvel = (pvapre.vel + pvacur.vel) / 2;
	qne_pre = Earth::qne(pvapre.pos);
	qne_cur = Earth::qne(pvacur.pos);
	temp1 = Rotation::quaternion2vector(qne_cur.inverse() * qne_pre);      // approximately q_downn(m-1)_upn(m)
	qne_mid = qne_pre * Rotation::rotvec2quaternion(temp1 / 2).inverse();  
	midpos[2] = (pvacur.pos[2] + pvapre.pos[2]) / 2;
	midpos = Earth::blh(qne_mid, midpos[2]);

	// 2. recompute rmrn, wie_n, wen_n at m-1/2
    Vector2d rmrn;
    Vector3d wie_n, wen_n;
	rmrn = Earth::meridianPrimeVerticalRadius(midpos[0]);
	wie_n << WGS84_WIE * cos(midpos[0]), 0, -WGS84_WIE * sin(midpos[0]);
	wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
		-midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

	// 3. n-frame rotation vector (n(m-1) with respect to n(m)-frame)
	temp1 = -(wie_n + wen_n) * ts;              // Gongmin Yan (4.1.11) from time m to time m-1, a '-' is needed  
	qnn = Rotation::rotvec2quaternion(temp1);   // q_downn(m-1)_upn(m)

	// 4. b-frame rotation vector (b(m) with respect to b(m-1)-frame)
//	temp1 = imucur.dtheta + imupre.dtheta.cross(imucur.dtheta) / 12;  
	temp1 = imucur.segment<3>(1) + imupre.segment<3>(1).cross(imucur.segment<3>(1)) / 12;  // Gongmin Yan (2.5.37)
	qbb = Rotation::rotvec2quaternion(temp1);

	// 5. attitude update finish         // Gongmin Yan (4.1.8) in quaternion form
	pvacur.qbn = qnn * pvapre.qbn * qbb; // q_down(m)_upb(m): q_downn(m-1)_upn(m) * q_downb(m-1)_upn(m-1) * q_downb(m)_upb(m-1)
	pvacur.Cbn = Rotation::quaternion2matrix(pvacur.qbn);   // DCM
	pvacur.att = Rotation::matrix2euler(pvacur.Cbn);        // Euler Angles
}