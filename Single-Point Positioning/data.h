#pragma once
#define c 2.99792458e8
#define mu_GPS 3.986005e14
#define mu_BDS 3.986004418e14
#define OMEGA_GPS 7.2921151467e-5
#define OMEGA_BDS 7.2921150e-5
#define Max 99999
#define X_WUHN -2267750.0744
#define Y_WUHN 5009154.1821
#define Z_WUHN 3221290.5515
#define R 6378000//m
#define h 375000//m
#define T0 20//C
#define P0 1013.25//mbar
#define RH0 0.5

#include<iostream>
#include <iomanip>
#include<math.h>
#include<string>
#include <fstream>
#include "coordinates.h"
#include"time.h"
#include<vector>
using namespace std;
//标准化数据
double DataNormalize(double input) {
	if (input == 0) return 0; 
	double absInput = fabs(input);
	
	int scale = static_cast<int>(log10(absInput));

	double normalizedData = absInput / pow(10, scale - 6);

	return (input < 0) ? -normalizedData : normalizedData;
}
double my_stod(const std::string& s) {
	for (char C : s) {
		if (!std::isspace(C)) {
			return std::stod(s);
		}
	}
	return 0.0;
}
class IONOCORR {
public:
	//satellite transmitted terms
	vector<double>TEMP;
	double GPSAlpha_n[4];//4系数
	double GPSBeta_n[4];
	double BDSAlpha_n[4];//4系数
	double BDSBeta_n[4];
};
class GMNREC {
public:
	string PRN;
	COMMONTIME TOC;
	double ClkBias;//a0
	double ClkDrift;//a1
	double ClkDriftRate;//a2

	double IODE;//数据龄期
	double Crs;//轨道半径的正弦调和改正项的振幅(m)
	double DeltaN;//卫星平均运动速率与计算值之差(PI/s)
	double M0;//参考时间的平近点角(PI)

	double Cuc;//纬度幅角的余弦调和改正项的振幅(rad)
	double e;//偏心率
	double Cus;//纬度幅角的正弦调和改正项的振幅(rad)
	double SqrtA;//长半轴的平方根(m^1/2)

	double TOE;//星历参考时间 (s)
	double Cic;//轨道倾角的余弦调和改正项的振幅(rad)
	double Omega;//按参考时间计算的升交点经度(PI)
	double Cis;//轨道倾角的正弦调和改正项的振幅(rad)

	double i0;//参考时间的轨道倾角(PI)
	double Crc;//轨道半径的余弦调和改正项的振幅(m)
	double omega;//近地点幅角(PI)
	double OmegaDot;//升交点赤经变化率(PI/s)

	double iDot;//轨道倾角变化率(PI/𝑠 )
	double CodesOnL2Channel;//CA码，以上数据足以计算单点定位
	double GPSWeek;
	double L2PDataFlag;//L2P码

	double SVAccuracy;//卫星精度
	double SVHealth;//健康值
	double TGD;//电离层延迟
	double IODC;//数据质量

	double TransTimeOfMsg;
	double Spare1;//拟合区间

	GMNREC() : PRN(""), ClkBias(0), ClkDrift(0), ClkDriftRate(0), IODE(0), Crs(0), DeltaN(0), M0(0),
		Cuc(0), e(0), Cus(0), SqrtA(0), TOE(0), Cic(0), Omega(0), Cis(0), i0(0), Crc(0),
		omega(0), OmegaDot(0), iDot(0), CodesOnL2Channel(0), GPSWeek(0), L2PDataFlag(0),
		SVAccuracy(0), SVHealth(0), TGD(0), IODC(0), TransTimeOfMsg(0), Spare1(0) {}
	bool isGPS() {
		if (PRN[0] == 'G')return 1;
		else return 0;
	}
	bool isBDS() {
		if (PRN[0] == 'C')return 1;
		else return 0;
	}
	bool isGEO() {
		if (PRN[0] == 'C' && ((stoi(PRN.substr(1, 2)) >= 1 && stoi(PRN.substr(1, 2)) <= 5) || (stoi(PRN.substr(1, 2)) >= 59 && stoi(PRN.substr(1, 2)) <= 60)))
			return 1;
		else return 0;
	}
};
class SATELLITE {
public:
	GMNREC gmn;
	
	double E;//高度角
	double A;//方位角
	
	double lambda;//经度
	double phi;//纬度
	//坐标
	CRDCARTESIAN xyz;

	//速度（地固坐标系）
	CRDCARTESIAN xyzDot;
	
	double DEL_tsv;//钟差
	double DEL_cld;//钟漂
	double DEL_tsvDot;	//钟速
	//电离层延迟
	double T_IONO;
	//对流层延迟
	double DEL_trop;
	//伪距
	double rho;
	//𝑡时刻测站与卫星之间的几何距离
	double RO;
	//卫星高度角和方位角计算
	void GetEA(CRDCARTESIAN rs, CRDCARTESIAN rr) {//卫星位置 rs(xyz)，接收机位置rr(xyz)
		CRDGEODETIC rsblh; 
		Cartesian2Geodetic(rs, rsblh);
		CRDGEODETIC rrblh; 
		Cartesian2Geodetic(rr, rrblh);
		double d = sqrt(pow(rs.x - rr.x, 2) + pow(rs.y - rr.y, 2) + pow(rs.z - rr.z, 2));
		double B = rrblh.longitude * PI / 180, L = rrblh.latitude * PI / 180;

		double h_[9] = { -sin(L),cos(L),0,
		-sin(B) * cos(L),-sin(B) * sin(L),cos(B),
		cos(B) * cos(L),cos(B) * sin(L),sin(B) };

		double e[3] = { rs.x - rr.x,rs.y - rr.y,rs.z = rr.z };
		Matrix H(3, 3, h_);
		Matrix E_(3, 1, e);
		Matrix ER;
		ER = H * E_;
		E = atan2(ER.data[0][0], ER.data[1][0]);
		A = asin(ER.data[2][0] / d);
	}
	void GetOrbNClk(COMMONTIME t) {//t为信号发射时刻的时间
		if (gmn.isGPS()) {
			//坐标计算
			double n0 = sqrt(mu_GPS / pow(gmn.SqrtA, 3));
			double tk = t.COM2GPSTime() - gmn.TOE;
			tk = (tk > 302400) ? (tk - 604800) : ((tk < -302400) ? (tk + 604800) : tk);
			double n = n0 + gmn.DeltaN;
			double Mk = gmn.M0 + n * tk;
			double Ek = Mk;
			double temp = 0;
			int iTimer = 0;
			do {
				temp = Ek;
				Ek = Mk + gmn.e * sin(Ek);
				iTimer++;
			} while ((fabs(Ek - temp) > 1e-12) &&( iTimer < 200));
			Mk = Ek - gmn.e * sin(Ek);
			double vk = atan2(sqrt(1 - gmn.e * gmn.e) * sin(Ek), cos(Ek) - gmn.e);
			double PHI_k = vk + gmn.omega;
			double DEL_uk = gmn.Cus * sin(2.0 * PHI_k) + gmn.Cuc * cos(2.0 * PHI_k);//δuk
			double DEL_rk = gmn.Crs * sin(2.0 * PHI_k) + gmn.Crc * cos(2.0 * PHI_k);
			double DEL_ik = gmn.Cis * sin(2.0 * PHI_k) + gmn.Cic * cos(2.0 * PHI_k);
			double uk = PHI_k + DEL_uk;
			double rk = gmn.SqrtA * gmn.SqrtA * (1 - gmn.e * cos(Ek)) + DEL_rk;
			double ik = gmn.i0 + DEL_ik + gmn.iDot * tk;
			double xk = rk * cos(uk);
			double yk = rk * sin(uk);
			//14-15
			double OmegaK = gmn.Omega + (gmn.OmegaDot - OMEGA_GPS) * tk - OMEGA_GPS * gmn.TOE;
			xyz.x = xk * cos(OmegaK) - yk * cos(ik) * sin(OmegaK);
			xyz.y = xk * sin(OmegaK) + yk * cos(ik) * cos(OmegaK);
			xyz.z = yk * sin(ik);

			////卫星速度计算
			double EkDot = n / (1 - gmn.e * cos(Ek));
			double PHI_kDot = sqrt((1 + gmn.e) / (1 - gmn.e)) * (cos(vk) / cos(Ek)) * (cos(vk) / cos(Ek)) * EkDot;
			double ukDot = PHI_kDot * (1 + 2 * gmn.Cus * cos(2 * PHI_k) - 2 * gmn.Cuc * sin(2 * PHI_k));
			double rkDot = EkDot * gmn.SqrtA * gmn.SqrtA * gmn.e * sin(Ek) + 2 * PHI_k * (gmn.Crs * cos(2 * PHI_k) - gmn.Crc * sin(2 * PHI_k));
			double ikDot = gmn.iDot + 2 * PHI_kDot * (gmn.Cis * cos(2 * PHI_k) - gmn.Cic * sin(2 * PHI_k));
			double xkDot = rkDot * cos(uk) - rk * ukDot * sin(uk);
			double ykDot = rkDot * sin(uk) + rk * ukDot * cos(uk);
			//5-6
			double OmegakDot = gmn.OmegaDot - OMEGA_GPS;
			xyzDot.x = xkDot * cos(PHI_k) - ykDot * cos(ik) * sin(PHI_k) - PHI_kDot * (xk * sin(PHI_k) + yk * cos(ik) * cos(PHI_k)) + ikDot * yk * sin(ik) * sin(PHI_k);
			xyzDot.y = xkDot * sin(PHI_k) + ykDot * cos(ik) * cos(PHI_k) + PHI_kDot * (xk * cos(PHI_k) - yk * cos(ik) * sin(PHI_k)) - ikDot * yk * sin(ik) * cos(PHI_k);
			xyzDot.z = ykDot * sin(ik) + ikDot * yk * cos(ik);

			//卫星钟差计算
			double DEL_tr = -2 * sqrt(mu_GPS) / (c * c) * gmn.e * gmn.SqrtA * sin(Ek);
			DEL_tsv = gmn.ClkBias + gmn.ClkDrift * (t.COM2GPSTime() - gmn.TOC.COM2GPSTime()) + gmn.ClkDriftRate * (t.COM2GPSTime() - gmn.TOC.COM2GPSTime()) * (t.COM2GPSTime() - gmn.TOC.COM2GPSTime()) + DEL_tr;
			//卫星钟漂计算
			DEL_cld = gmn.ClkDrift + 2 * gmn.ClkDriftRate * fabs(gmn.TOC - t);
			//卫星钟速计算
			double DEL_trDot = -2 * sqrt(mu_GPS) / (c * c) * gmn.e * gmn.SqrtA * cos(Ek) * EkDot;
			DEL_tsvDot = gmn.ClkDrift + 2 * gmn.ClkDriftRate * (t.COM2GPSTime() - gmn.TOC.COM2GPSTime()) + DEL_trDot;
		}
		else {
			//坐标计算
			double n0 = sqrt(mu_GPS / pow(gmn.SqrtA, 3));
			double tk = t.COM2GPSTime() - gmn.TOE;
			tk = (tk > 302400) ? (tk - 604800) : ((tk < -302400) ? (tk + 604800) : tk);
			double n = n0 + gmn.DeltaN;
			double Mk = gmn.M0 + n * tk;
			double Ek = Mk;
			double temp = 0;
			int iTimer = 0;
			do {
				temp = Ek;
				Ek = Mk + gmn.e * sin(Ek);
				iTimer++;
			} while (fabs(Ek - temp) > 1e-12 && iTimer < 200);
			Mk = Ek - gmn.e * sin(Ek);
			double vk = atan2(sqrt(1 - gmn.e * gmn.e) * sin(Ek), cos(Ek) - gmn.e);
			double PHI_k = vk + gmn.omega;
			double DEL_uk = gmn.Cus * sin(2.0 * PHI_k) + gmn.Cuc * cos(2.0 * PHI_k);//δuk
			double DEL_rk = gmn.Crs * sin(2.0 * PHI_k) + gmn.Crc * cos(2.0 * PHI_k);
			double DEL_ik = gmn.Cis * sin(2.0 * PHI_k) + gmn.Cic * cos(2.0 * PHI_k);
			double uk = PHI_k + DEL_uk;
			double rk = gmn.SqrtA * gmn.SqrtA * (1 - gmn.e * cos(Ek)) + DEL_rk;
			double ik = gmn.i0 + DEL_ik + gmn.iDot * tk;
			double xk = rk * cos(uk);
			double yk = rk * sin(uk);
			double OmegaK, XGK, YGK, ZGK;

			if (gmn.isGEO()) {
				OmegaK = gmn.Omega + gmn.OmegaDot * tk - OMEGA_GPS * gmn.TOE;
				XGK = xk * cos(OmegaK) - yk * cos(ik) * sin(OmegaK);
				YGK = xk * sin(OmegaK) + yk * cos(ik) * cos(OmegaK);
				ZGK = yk * sin(ik);
				xyz.x = cos(OMEGA_BDS * tk) * XGK + sin(OMEGA_BDS * tk) * cos(-5 / 180 * PI) * YGK + sin(OMEGA_BDS * tk) * sin(-5 / 180 * PI) * ZGK;
				xyz.y = -sin(OMEGA_BDS * tk) * XGK + cos(OMEGA_BDS * tk) * cos(-5 / 180 * PI) * YGK + cos(OMEGA_BDS * tk) * sin(-5 / 180 * PI) * ZGK;
				xyz.z = -sin(-5 / 180 * PI) * YGK + cos(-5 / 180 * PI) * ZGK;
			}
			else {
				OmegaK = gmn.Omega + (gmn.OmegaDot - OMEGA_GPS) * tk - OMEGA_GPS * gmn.TOE;
				xyz.x = xk * cos(OmegaK) - yk * cos(ik) * sin(OmegaK);
				xyz.y = xk * sin(OmegaK) + yk * cos(ik) * cos(OmegaK);
				xyz.z = yk * sin(ik);
			}

			//卫星速度计算
			double EkDot = n / (1 - gmn.e * cos(Ek));
			double PHI_kDot = sqrt((1 + gmn.e) / (1 - gmn.e)) * (cos(vk) / cos(Ek)) * (cos(vk) / cos(Ek)) * EkDot;
			double ukDot = PHI_kDot * (1 + 2 * gmn.Cus * cos(2 * PHI_k) - 2 * gmn.Cuc * sin(2 * PHI_k));
			double rkDot = EkDot * gmn.SqrtA * gmn.SqrtA * gmn.e * sin(Ek) + 2 * PHI_k * (gmn.Crs * cos(2 * PHI_k) - gmn.Crc * sin(2 * PHI_k));
			double ikDot = gmn.iDot + 2 * PHI_kDot * (gmn.Cis * cos(2 * PHI_k) - gmn.Cic * sin(2 * PHI_k));
			double xkDot = rkDot * cos(uk) - rk * ukDot * sin(uk);
			double ykDot = rkDot * sin(uk) + rk * ukDot * cos(uk);
			double OmegakDot;
			//5-6
			if (gmn.isGEO()) {
				OmegakDot = gmn.OmegaDot;
				double XGKDot = xkDot * cos(PHI_k) - ykDot * cos(ik) * sin(PHI_k) - PHI_kDot * (xk * sin(PHI_k) + yk * cos(ik) * cos(PHI_k)) + ikDot * yk * sin(ik) * sin(PHI_k);
				double YGKDot = xkDot * sin(PHI_k) + ykDot * cos(ik) * cos(PHI_k) + PHI_kDot * (xk * cos(PHI_k) - yk * cos(ik) * sin(PHI_k)) - ikDot * yk * sin(ik) * cos(PHI_k);
				double ZGKDot = ykDot * sin(ik) + ikDot * yk * cos(ik);
				xyzDot.x = OMEGA_BDS * (-sin(OMEGA_BDS * tk) * XGK + cos(OMEGA_BDS * tk) * cos(-5 / 180 * PI) * YGK + cos(OMEGA_BDS * tk) * sin(-5 / 180 * PI) * ZGK) + cos(OMEGA_BDS * tk) * XGKDot + sin(OMEGA_BDS * tk) * cos(-5 / 180 * PI) * YGKDot + sin(OMEGA_BDS * tk) * sin(-5 / 180 * PI) * ZGKDot;
				xyzDot.y = OMEGA_BDS * (-cos(OMEGA_BDS * tk) * XGK - sin(OMEGA_BDS * tk) * cos(-5 / 180 * PI) * YGK - sin(OMEGA_BDS * tk) * sin(-5 / 180 * PI) * ZGK) - sin(OMEGA_BDS * tk) * XGKDot + cos(OMEGA_BDS * tk) * cos(-5 / 180 * PI) * YGKDot + cos(OMEGA_BDS * tk) * sin(-5 / 180 * PI) * ZGKDot;
				xyzDot.z = -sin(-5 / 180 * PI) * YGKDot + cos(-5 / 180 * PI) * ZGKDot;
			}
			else {
				OmegakDot = gmn.OmegaDot - OMEGA_GPS;
				xyzDot.x = xkDot * cos(PHI_k) - ykDot * cos(ik) * sin(PHI_k) - PHI_kDot * (xk * sin(PHI_k) + yk * cos(ik) * cos(PHI_k)) + ikDot * yk * sin(ik) * sin(PHI_k);
				xyzDot.y = xkDot * sin(PHI_k) + ykDot * cos(ik) * cos(PHI_k) + PHI_kDot * (xk * cos(PHI_k) - yk * cos(ik) * sin(PHI_k)) - ikDot * yk * sin(ik) * cos(PHI_k);
				xyzDot.z = ykDot * sin(ik) + ikDot * yk * cos(ik);
			}

			//卫星钟差计算
			double DEL_tr = -2 * sqrt(mu_BDS) / (c * c) * gmn.e * gmn.SqrtA * sin(Ek);
			DEL_tsv = gmn.ClkBias + gmn.ClkDrift * (t.COM2GPSTime() - gmn.TOC.COM2GPSTime()) + gmn.ClkDriftRate * (t.COM2GPSTime() - gmn.TOC.COM2GPSTime()) * (t.COM2GPSTime() - gmn.TOC.COM2GPSTime()) + DEL_tr;
			//卫星钟速计算
			double DEL_trDot = -2 * sqrt(mu_BDS) / (c * c) * gmn.e * gmn.SqrtA * cos(Ek) * EkDot;
			DEL_tsvDot = gmn.ClkDrift + 2 * gmn.ClkDriftRate * (t.COM2GPSTime() - gmn.TOC.COM2GPSTime()) + DEL_trDot;


		}
		// 计算站点的经度和纬度
		lambda = atan2(Y_WUHN, X_WUHN); // 经度
		phi = atan2(Z_WUHN, sqrt(X_WUHN * X_WUHN + Y_WUHN * Y_WUHN)); // 纬度
		//计算伪距
		rho = fabs(gmn.TOC - t) * c;
	/*	cout << rho << "  ";*/
		int l = 0;
		
	
	}
	void PrintOrbNClk() {
		cout << "PRN:" << gmn.PRN << endl;
		cout << "坐标:" << endl;
		cout << "Xk:" << fixed << setprecision(9) << xyz.x << endl;
		cout << "Yk:" << fixed << setprecision(9) << xyz.y << endl;
		cout << "Zk:" << fixed << setprecision(9) << xyz.z << endl;
		//cout << "速度:" << endl;
		//cout << "XkDot:" << fixed << setprecision(9) << xyzDot.x << endl;
		//cout << "YkDot:" << fixed << setprecision(9) << xyzDot.y << endl;
		//cout << "ZkDot:" << fixed << setprecision(9) << xyzDot.z << endl;
		cout << "钟差:" << endl;
		cout << "DEL_tsv:" << fixed << setprecision(9) << DEL_tsv << endl;
		/*cout << "DEL_tsvDot:" << fixed << setprecision(9) << DEL_tsvDot << endl;*/
	}
	//电离层延迟
	void GetIONO(COMMONTIME t_,IONOCORR corr) {
		GPSTime gpst;
		GPSTime2CommonTime(t_, gpst);
		if (gmn.isGPS()) {
			double PHI = 0.0137 / (E + 0.11) - 0.022;//semi-circles

			double PHI_i = phi + PHI * cos(A);
			if (PHI_i > 0.416)PHI_i = 0.416;
			else if (PHI_i < -0.416)PHI_i = -0.416;

			double MU_i = lambda + PHI * sin(A) / cos(PHI_i);

			double PHI_m = PHI_i + 0.064 * cos(MU_i - 1.617);

			double t = 4.32e4 * PHI_i + gpst.sn;
			if (t > 86400)t -= 86400;
			else if (t < 0)t += 86400;

			double F = 1 + 16.0 * pow((0.53 - E), 3);

			double PER = 0;
			for (int i = 0; i < 4; i++) {
				PER += corr.GPSBeta_n[i] * pow(PHI_m, i);
			}
			if (PER < 72000)PER = 72000;//sec

			double AMP = 0;
			for (int i = 0; i < 4; i++) {
				AMP += corr.GPSAlpha_n[i] * pow(PHI_m, i);
			}
			if (AMP < 0)AMP = 0;

			double x = 2 * PI * (t - 50400) / PER;//radians

			T_IONO = (fabs(x) < 1.57) ? F * (5.0e-9 + AMP * (1 - x * x / 2 + x * x * x * x / 24)) : F * 5.0e-9;//sec
		}
		else if (gmn.isBDS()) {
			double PHI = PI / 2 - E - asin(R * cos(E) / (R + h));

			double PHI_M = asin(sin(phi) * cos(PHI) + cos(phi) * sin(PHI) * cos(A));

			double MU_M = lambda + asin(sin(PHI) * sin(A) / cos(PHI_M));

			double A2 = 0;
			for (int i = 0; i < 4; i++) {
				A2 += corr.BDSAlpha_n[i] * pow(PHI_M / PI, i);
			}
			if (A2 < 0)A2 = 0;

			double A4 = 0;
			for (int i = 0; i < 4; i++) {
				A4 += corr.BDSBeta_n[i] * pow(PHI_M / PI, i);
			}
			if (A2 >= 172800)A2 = 172800;
			else if (A2 < 72000)A4 = 72000;

			double t = fmod((gpst.sn + MU_M * 43200 / PI), 86400);

			double T_IONO_z = 5e-9 + A2 * cos(2 * PI * (t - 50400) / A4);
			if (fabs(t - 50400) >= A4 / 4)T_IONO_z = 5e-9;

			T_IONO = T_IONO_z / sqrt(1 - (R * cos(E) / (R + h)) * (R * cos(E) / (R + h)));
		}
		
	}
	//对流层延迟 h的范围要求：适用于海拔在0到约11公里之间。h对37594取模

	void GetTROP(CRDCARTESIAN site_) {
		CRDGEODETIC site;
		Cartesian2Geodetic(site_, site);
		double T = T0 - 0.065 * site.height;
		double P = P0 * pow((1 - 0.0000266 * fmod(site.height,37594)), 5.225);
		double RH = pow(RH0, -0.0006369 * site.height);
		double e = pow(RH, -37.2465 + 0.213166 * T - 0.000256908 * T);

		double hd = 40136 + 148.72 * (T - 273.16);
		double hw = 11000.0;
		double Kw = 155.2e-7 * 4810 * e * (hw - site.height) / (T * T);
		double Kd = 155.2e-7 * P * (hw - site.height) / T;
		DEL_trop = Kd / sin(sqrt(E * E + 6.25)) + Kw / sin(sqrt(E * E + 2.25));
	}
};

