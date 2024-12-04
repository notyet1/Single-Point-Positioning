#pragma once
#include<iostream>
#include<string>
#include <algorithm> 
#include<math.h>
#include<ctype.h>
#include <sstream>
#include <fstream>
#include"time.h"
#include"matrix.h"
#include"coordinates.h"
#include"data.h"
#define OMEGA 7.2921151467e-5
#define c 2.99792458e8
class DATA {
public:
	COMMONTIME ts;//信号发射时刻
	COMMONTIME tr;//信号接收时刻
	SATELLITE *sat;
	int num;
	int nums;
	int gpsnum;
	int bdsnum;
	CRDCARTESIAN Station;//测站坐标 APPROX POSITION XYZ 
	CRDCARTESIAN StationSpeed;//测站速度
	//接收机钟差δ𝑡𝑖(m)
	double bu_gps;
	double bu_bds;
	//接收机钟漂δ𝑡𝑖Dot(m/s)
	double bd_gps;
	double bd_bds;
	IONOCORR corr;
	//最小二乘 V=B*X-W
	vector<double>l_, m_, n_,B_,W_,P_;
	Matrix l, m, n;//方向余弦S*1
	Matrix W;//OMC向量矩阵S*1
	Matrix B;//参数系数矩阵S*4
	Matrix P;//权重经验矩阵P=sinE S*S
	Matrix Q;//计算矩阵4*4
	Matrix Q_;
	Matrix X;//4*1
	//精度因子
	double GDOP;
	double PDOP;
	//已选择的卫星数组
	vector<SATELLITE>sel;
	vector<SATELLITE>GPSS;
	vector<SATELLITE>BDSS;
	DATA():num(2880),tr(2024,7,27,0,0,0),Station(-2267796.9641, 5009421.6975, 3220952.5436) ,StationSpeed(0,0,0),bu_gps(0), bu_bds(0),gpsnum(0),bdsnum(0) {
		sat = new SATELLITE [num];
	}
	~DATA() {
		delete [] sat;
	}
	//选择卫星
	void SelectSatellite(string s) {
		std::istringstream stream(s);
		int number;
		while (stream >> number) {
			sel.push_back(sat[number-1]);
			nums++;
		}
	}
	//卫星分类计数
	void SatelliteCount(vector<SATELLITE>sel) {
		for (SATELLITE SAT : sel) {
			if (SAT.gmn.isGPS()) {
				gpsnum++;
			}
			else if (SAT.gmn.isBDS()) {
				bdsnum++;
			}
		}
	}
	//卫星分类
	void SatelliteClassification(vector<SATELLITE>sel) {
		for (SATELLITE SAT : sel) {
			if (SAT.gmn.isGPS()) {
				GPSS.push_back(SAT);
			}
			else if (SAT.gmn.isBDS()) {
				BDSS.push_back(SAT);
			}
		}
	}
	//确定定位方式
	bool isSinglePosition(vector<SATELLITE> sel) {
		int hasgps = false;
		int hasbds = false;
		for (SATELLITE sat : sel) {
			if (sat.gmn.isGPS())hasgps = true;
			else if(sat.gmn.isBDS())hasbds = true;
			if (hasgps && hasbds) return false;
		}
		return true;
	}
	//初始化
	void DATA_INIT_SSP() {
		l_.resize(nums); m_.resize(nums); n_.resize(nums);
		B_.resize(nums*4); W_.resize(nums); P_.resize(nums*nums);
	
		/*l_.clear(); m_.clear(); n_.clear(),B_.clear(),W_.clear(),P_.clear();
		W.clear(); B.clear(); P.clear(); Q.clear(); X.clear();*/
	}
	void DATA_INIT_CSP() {
		l_.resize(nums); m_.resize(nums); n_.resize(nums);
		B_.resize(nums * 5); W_.resize(nums); P_.resize(nums * nums);
	}

	//单系统单点定位
	void GetSSPP() {
		vector<double> xUpdates;
		vector<double> yUpdates;
		vector<double> zUpdates;
		ofstream outFile("XYZerror.csv");
		outFile << "X_error,Y_error,Z_error" << std::endl;
		if (nums < 4) {
			cout << "卫星数量太少，无法定位" << endl;
		}
		else if (nums >= 4) {
			if(sel[0].gmn.isGPS()){
				int Iteration = 0;
				while (1) {
					Iteration++;
					/*cout << Iteration << endl;*/
					
					double alpha;
					//迭代计算参数,地卫距𝜌与卫星位置有关，需要迭代计算
					for (int i = 0; i < nums; i++) {
						//信号传播时间初值
						double Del_t = 0.075;
						double Iteration_count = 0;
						
						while (1) {
							Iteration_count++;
							double Del_t_ = Del_t;
							double d = tr.totalDays() - (bu_gps / c) / 86400.0 - Del_t;
							ts = JulianDay2COM2024(d);
							
							sel[i].GetOrbNClk(ts);

							//地球自转改正
							alpha = OMEGA * sel[i].rho / c;
							
							sel[i].xyz.x = sel[i].xyz.x * cos(alpha) + sel[i].xyz.y * sin(alpha);
							sel[i].xyz.y = -sel[i].xyz.x * sin(alpha) + sel[i].xyz.y * cos(alpha);
							double D= (sel[i].xyz.x - Station.x) * (sel[i].xyz.x - Station.x) 
							+ (sel[i].xyz.y - Station.y) * (sel[i].xyz.y - Station.y);
							 +(sel[i].xyz.z - Station.z) * (sel[i].xyz.z - Station.z);
							sel[i].RO = sqrt(D);
							Del_t = sel[i].RO / c;
							bu_gps = sel[i].rho - sel[i].RO;//???
							if (fabs(Del_t - Del_t_) <= 1e-12 || Iteration_count >= 6) {
								
								break;
							}
						}
						sel[i].GetEA(sel[i].xyz, Station);
						sel[i].GetIONO(ts, corr);
						sel[i].GetTROP(Station);
						
					}

					//设置方向余弦矩阵l m n S*1
					for (int i = 0; i < nums; i++) {
						l_[i] = (sel[i].xyz.x - Station.x) / sel[i].RO;
						m_[i] = (sel[i].xyz.y - Station.y) / sel[i].RO;
						n_[i] = (sel[i].xyz.z - Station.z) / sel[i].RO;
					}                                                                                
					l.setMatrix(nums, 1, l_.data());
					m.setMatrix(nums, 1, m_.data());
					n.setMatrix(nums, 1, n_.data());
					//设置参数系数矩阵B S*4
					for (int i = 0; i < nums; i++) {
						B_[4 * i + 0] = -l_[i];
						B_[4 * i + 1] = -m_[i];
						B_[4 * i + 2] = -n_[i];
						B_[4 * i + 3] = 1;
					}
					B.setMatrix(nums, 4, B_.data());

					//if (Iteration == 2) {/////////////////
					//	int tip=1;
					//}

					//设置OMC向量矩阵W S*1
					for (int i = 0; i < nums; i++) {
						W_[i] = sel[i].rho - sel[i].RO - bu_gps + c * sel[i].DEL_tsv - sel[i].T_IONO - sel[i].DEL_trop;
					}
					W.setMatrix(nums, 1, W_.data());
					//设置权重经验矩阵P S*S
					for (int i = 0; i < nums; i++) {
						P_[i * nums + i] = sin(sel[i].E);
					}
					P.setMatrix(nums, nums, P_.data());
					//设置计算矩阵Q 4*4(用于精度评定)
					Q_ = B.transpose() * B;
					Q_ = Q_.regularizeMatrix(1e-8);
					Q = Q_.inverse();
					//最小二乘
					X = reweightedLeastSquares(B, W, P);
					//接收机坐标，钟差数值更新
					double x_temp= Station.x, y_temp= Station.y, z_temp= Station.z;

					Station.x += DataNormalize(X.data[0][0]);
					Station.y += DataNormalize(X.data[1][0]);
					Station.z += DataNormalize(X.data[2][0]);
					bu_gps += DataNormalize(X.data[3][0]);
					xUpdates.push_back(Station.x-x_temp);
					yUpdates.push_back(Station.y-y_temp);
					zUpdates.push_back(Station.z-z_temp);
					/*cout << "迭代次数：" << Iteration<< endl;
					cout << Station.x << " " << Station.y << " " << Station.z;*/
					//设置最大迭代次数7，精度评定
					if (((abs(X.data[0][0]) < 1e-4 && abs(X.data[1][0]) < 1e-4 && abs(X.data[2][0]) < 1e-4)) || Iteration >= 7) {
						GDOP = sqrt(Q.trace());
						PDOP = sqrt(Q.trace() - Q.data[3][3] * Q.data[3][3]);
						break;
					}
				}
				for (size_t i = 0; i < xUpdates.size(); ++i) {
					outFile << xUpdates[i] << "," << yUpdates[i] << "," << zUpdates[i] << std::endl;
				}

				outFile.close();
				//输出定位坐标、接收机钟差及其精度信息
				cout << "定位坐标:" << Station.x << "," << Station.y << "," << Station.z << endl;
				cout << "接收机钟差:" << bu_gps << endl;
				cout << "精度信息 GDOP:" << GDOP << " PDOP:" << PDOP << endl;
			}
			else if(sel[0].gmn.isBDS()) {
				int Iteration = 0;
				while (1) {
					Iteration++;
					double alpha;

					//迭代计算参数
					for (int i = 0; i < nums; i++) {

						double Del_t = 0.075;
						double Iteration_count = 0;
						while (1) {
							Iteration_count++;
							double Del_t_ = Del_t;
							ts = JulianDay2COM2024(tr.calculateJulianDay() - bu_bds / c - Del_t);
							sel[i].GetOrbNClk(ts);
							//地球自转改正
							alpha = OMEGA * sel[i].rho / c;
							sel[i].xyz.x = sel[i].xyz.x * cos(alpha) + sel[i].xyz.y * sin(alpha);
							sel[i].xyz.y = -sel[i].xyz.x * sin(alpha) + sel[i].xyz.y * cos(alpha);
							sel[i].RO = sqrt((sel[i].xyz.x - Station.x) * (sel[i].xyz.x - Station.x) + (sel[i].xyz.y - Station.y) * (sel[i].xyz.y - Station.y) + (sel[i].xyz.z - Station.z) * (sel[i].xyz.z - Station.z));
							Del_t = sel[i].RO / c;
							bu_bds = sel[i].rho - sel[i].RO;//???
							if (fabs(Del_t - Del_t_) <= 1e-12 || Iteration_count >= 6)break;
						}
						sel[i].GetEA(sel[i].xyz, Station);
						sel[i].GetIONO(ts, corr);
						sel[i].GetTROP(Station);
						
					}
					//设置方向余弦矩阵l m n S*1
					for (int i = 0; i < nums; i++) {
						l_[i] = (sel[i].xyz.x - Station.x) / sel[i].RO;
						m_[i] = (sel[i].xyz.y - Station.y) / sel[i].RO;
						n_[i] = (sel[i].xyz.z - Station.z) / sel[i].RO;
					}
					l.setMatrix(nums, 1, l_.data());
					m.setMatrix(nums, 1, m_.data());
					n.setMatrix(nums, 1, n_.data());
					//设置参数系数矩阵B S*4
					for (int i = 0; i < nums; i++) {
						B_[4 * i + 0] = -l_[i];
						B_[4 * i + 1] = -m_[i];
						B_[4 * i + 2] = -n_[i];
						B_[4 * i + 3] = 1;
					}
					B.setMatrix(nums, 4, B_.data());
					//设置OMC向量矩阵W S*1
					for (int i = 0; i < nums; i++) {
						W_[i] = sel[i].rho - sel[i].RO - bu_bds + c * sel[i].DEL_tsv - sel[i].T_IONO - sel[i].DEL_trop;
					}
					W.setMatrix(nums, 1, W_.data());
					//设置权重经验矩阵P S*S
					for (int i = 0; i < nums; i++) {
						P_[i * nums + i] = sin(sel[i].E);
					}
					P.setMatrix(nums, nums, P_.data());
					//设置计算矩阵Q 4*4
					Q_ = B.transpose() * B;
					Q_ = Q_.regularizeMatrix(1e-8);
					/*Q = Q_.inverse();*/
					//最小二乘
					X = reweightedLeastSquares(B, W, P);
					double x_temp = Station.x, y_temp = Station.y, z_temp = Station.z;
					//接收机坐标，钟差数值更新
					Station.x += DataNormalize(X.data[0][0]);
					Station.y += DataNormalize(X.data[1][0]);
					Station.z += DataNormalize(X.data[2][0]);
					bu_bds += DataNormalize(X.data[3][0]);
					xUpdates.push_back(Station.x - x_temp);
					yUpdates.push_back(Station.y - y_temp);
					zUpdates.push_back(Station.z - z_temp);
					//设置最大迭代次数7，精度评定
					if (((abs(X.data[0][0]) < 1e-4 && abs(X.data[1][0]) < 1e-4 && abs(X.data[2][0]) < 1e-4)) || Iteration >= 7) {
						GDOP = sqrt(Q.trace());
						PDOP = sqrt(Q.trace() - Q.data[3][3] * Q.data[3][3]);
						break;
					}
				}
				for (size_t i = 0; i < xUpdates.size(); ++i) {
					outFile << xUpdates[i] << "," << yUpdates[i] << "," << zUpdates[i] << std::endl;
				}

				outFile.close();
				//输出定位坐标、接收机钟差及其精度信息
				cout << "定位坐标:" << Station.x << "," << Station.y << "," << Station.z << endl;
				cout << "接收机钟差:" << bu_bds << endl;
				cout << "精度信息 GDOP:" << GDOP << " PDOP:" << PDOP << endl;
			}
		}
	}
	//单系统单点测速
	void GetSSPST() {
		W_.clear();  W.clear();
		//计算地卫距
		for (int i = 0; i < nums; i++) {
			sel[i].RO = sqrt((sel[i].xyz.x - Station.x) * (sel[i].xyz.x - Station.x) + (sel[i].xyz.y - Station.y) * (sel[i].xyz.y - Station.y) + (sel[i].xyz.z - Station.z) * (sel[i].xyz.z - Station.z));
		}

		//设置OMC向量矩阵W
		for (int i = 0; i < nums; i++) {
			double w = ((sel[i].xyz.x - Station.x) * (sel[i].xyzDot.x - StationSpeed.x) + (sel[i].xyz.y - Station.y) * (sel[i].xyzDot.y - StationSpeed.y) + (sel[i].xyz.z - Station.z) * (sel[i].xyzDot.z - StationSpeed.z)) / sel[i].RO + l_[i] * sel[i].xyzDot.x + m_[i] * sel[i].xyzDot.y + n_[i] * sel[i].xyzDot.z + c * sel[i].DEL_cld;
			W_.push_back(w);
		}
		W.setMatrix(nums, 1, W_.data());
		//最小二乘
		X = reweightedLeastSquares(B, W, P);
		//接收机速度,钟漂数值更新,无须迭代计算
		StationSpeed.x += X.data[0][0];
		StationSpeed.y += X.data[1][0];
		StationSpeed.z += X.data[2][0];
		if (sel[0].gmn.isGPS()){
			bd_gps+= X.data[3][0];
			//精度评定
			GDOP = sqrt(Q.trace());
			PDOP = sqrt(Q.trace() - Q.data[3][3] * Q.data[3][3]);
			//输出接收机速度、钟漂及其精度信息
			cout << "接收机速度:" << StationSpeed.x << "," << StationSpeed.y << "," << StationSpeed.z << endl;
			cout << "接收机钟漂:" << bd_gps << endl;
			cout << "精度信息 GDOP:" << GDOP << " PDOP:" << PDOP << endl;
		}
		else  if (sel[0].gmn.isBDS()) {
			bd_bds += X.data[3][0];
			//精度评定
			GDOP = sqrt(Q.trace());
			PDOP = sqrt(Q.trace() - Q.data[3][3] * Q.data[3][3]);
			//输出接收机速度、钟漂及其精度信息
			cout << "接收机速度:" << StationSpeed.x << "," << StationSpeed.y << "," << StationSpeed.z << endl;
			cout << "接收机钟漂:" << bd_bds << endl;
			cout << "精度信息 GDOP:" << GDOP << " PDOP:" << PDOP << endl;
		}

	}
	//组合单点定位
	void GetCSPP() {		
		vector<double> xUpdates;
		vector<double> yUpdates;
		vector<double> zUpdates;
		ofstream outFile("XYZerror.csv");
		outFile << "X_error,Y_error,Z_error" << std::endl;
		if (nums < 4) {
			cout << "卫星数量太少，无法定位" << endl;
		}
		else if (nums >= 4) {
			int Iteration = 0;
			while (1) {
				Iteration++;
				
				double alpha;

				//迭代计算参数
				for (int i = 0; i < nums; i++) {
					
					double Del_t = 0.075;
					double Iteration_count = 0;
					if(sel[i].gmn.isGPS()){
						while (1) {
							Iteration_count++;
							
							double Del_t_ = Del_t;
							double d = (tr.totalDays()*86400 - bu_gps / c  - Del_t)/86400;
							ts = JulianDay2COM2024(d);
							sel[i].GetOrbNClk(ts);
							//地球自转改正
							alpha = OMEGA * sel[i].rho / c;
							sel[i].xyz.x = sel[i].xyz.x * cos(alpha) + sel[i].xyz.y * sin(alpha);
							sel[i].xyz.y = -sel[i].xyz.x * sin(alpha) + sel[i].xyz.y * cos(alpha);
							sel[i].RO = sqrt((sel[i].xyz.x - Station.x) * (sel[i].xyz.x - Station.x) + (sel[i].xyz.y - Station.y) * (sel[i].xyz.y - Station.y) + (sel[i].xyz.z - Station.z) * (sel[i].xyz.z - Station.z));
							Del_t = sel[i].RO / c;
							bu_gps = sel[i].rho - sel[i].RO;//???
							if (fabs(Del_t - Del_t_) <= 1e-12 || Iteration_count >= 6) { 
								
								break; 
							}
						}
					}
					else if (sel[i].gmn.isBDS()) {
						while (1) {
							Iteration_count++;
							double Del_t_ = Del_t;
							ts = JulianDay2COM2024(tr.calculateJulianDay() - bu_bds / c - Del_t);
							sel[i].GetOrbNClk(ts);
							//地球自转改正
							alpha = OMEGA * sel[i].rho / c;
							sel[i].xyz.x = sel[i].xyz.x * cos(alpha) + sel[i].xyz.y * sin(alpha);
							sel[i].xyz.y = -sel[i].xyz.x * sin(alpha) + sel[i].xyz.y * cos(alpha);
							sel[i].RO = sqrt((sel[i].xyz.x - Station.x) * (sel[i].xyz.x - Station.x) + (sel[i].xyz.y - Station.y) * (sel[i].xyz.y - Station.y) + (sel[i].xyz.z - Station.z) * (sel[i].xyz.z - Station.z));
							Del_t = sel[i].RO / c;
							bu_bds = sel[i].rho - sel[i].RO;//???
							if (fabs(Del_t - Del_t_) <= 1e-12 || Iteration_count >= 6)break;
						}
					}
					
					sel[i].GetEA(sel[i].xyz, Station);
					sel[i].GetIONO(ts, corr);
					sel[i].GetTROP(Station);
				
				}
				//更新GPSS和BDSS
				SatelliteClassification(sel);
				//设置方向余弦矩阵l m n
				for (int i = 0; i < gpsnum; i++) {
					l_[i] = (GPSS[i].xyz.x - Station.x) / GPSS[i].RO;
					m_[i] = (GPSS[i].xyz.y - Station.y) / GPSS[i].RO;
					n_[i] = (GPSS[i].xyz.z - Station.z) / GPSS[i].RO;
				}
				for (int i = gpsnum; i < nums; i++) {
					l_[i] = (BDSS[i - gpsnum].xyz.x - Station.x) / BDSS[i - gpsnum].RO;
					m_[i] = (BDSS[i - gpsnum].xyz.y - Station.y) / BDSS[i - gpsnum].RO;
					n_[i] = (BDSS[i - gpsnum].xyz.z - Station.z) / BDSS[i - gpsnum].RO;
				}
				l.setMatrix(nums, 1, l_.data());
				m.setMatrix(nums, 1, m_.data());
				n.setMatrix(nums, 1, n_.data());
				//设置参数系数矩阵B S*5
				for (int i = 0; i < nums; i++) {
					B_[4 * i + 0] = -l_[i];
					B_[4 * i + 1] = -m_[i];
					B_[4 * i + 2] = -n_[i];
				}
				for (int i = 0; i < gpsnum; i++) {
					B_[4 * i + 3] = 1; 
					B_[4 * i + 4] = 0;
				}
				for (int i = gpsnum; i < nums; i++) {
					B_[4 * i + 3] = 0;
					B_[4 * i + 4] = 1;
				}
				B.setMatrix(nums, 5, B_.data());
				//设置OMC向量矩阵W S*1
				for (int i = 0; i < gpsnum; i++) {
					W_[i] = GPSS[i].rho - GPSS[i].RO - bu_gps + c * GPSS[i].DEL_tsv - GPSS[i].T_IONO - GPSS[i].DEL_trop;
				}
				for (int i = gpsnum; i < nums; i++) {
					W_[i] = BDSS[i - gpsnum].rho - BDSS[i - gpsnum].RO - bu_bds + c * BDSS[i - gpsnum].DEL_tsv - BDSS[i - gpsnum].T_IONO - BDSS[i - gpsnum].DEL_trop;
				}
				W.setMatrix(nums, 1, W_.data());
				//设置权重经验矩阵P S*S
				for (int i = 0; i < gpsnum; i++) {
					P_[i * nums + i] = sin(GPSS[i].E);
				}
				for (int i = gpsnum; i < nums; i++) {
					P_[i * nums + i] = sin(BDSS[i - gpsnum].E);
				}
				P.setMatrix(nums, nums, P_.data());
				//设置计算矩阵Q
				Q = (B.transpose()*P * B).inverse();
				//最小二乘
				X = reweightedLeastSquares(B, W, P);
				double x_temp = Station.x, y_temp = Station.y, z_temp = Station.z;
				//接收机坐标，钟差数值更新
				Station.x += DataNormalize(X.data[0][0]);
				Station.y += DataNormalize(X.data[1][0]);
				Station.z += DataNormalize(X.data[2][0]);
				bu_gps += DataNormalize(X.data[3][0]);
				bu_bds += DataNormalize(X.data[4][0]);
				xUpdates.push_back(Station.x - x_temp);
				yUpdates.push_back(Station.y - y_temp);
				zUpdates.push_back(Station.z - z_temp);
				//设置最大迭代次数7，精度评定
				if (((abs(X.data[0][0]) < 1e-4 && abs(X.data[1][0]) < 1e-4 && abs(X.data[2][0]) < 1e-4)) || Iteration >= 7) {
					GDOP = sqrt(Q.trace());
					PDOP = sqrt(Q.trace() - Q.data[3][3] * Q.data[3][3]);
					break;
				}
			}
			for (size_t i = 0; i < xUpdates.size(); ++i) {
				outFile << xUpdates[i] << "," << yUpdates[i] << "," << zUpdates[i] << std::endl;
			}

			outFile.close();
			//输出定位坐标、接收机钟差及其精度信息
			cout << "定位坐标:" << Station.x << "," << Station.y << "," << Station.z << endl;
			cout << "接收机钟差 " <<"GPS系统:" << bu_gps << " 北斗系统:" <<bu_bds<<endl;
			cout << "精度信息 GDOP:" << GDOP << " PDOP:" << PDOP << endl;
		}
	}
	//组合单点测速
	void GetCSPST() {
		W_.clear(), W.clear();
		//计算地卫距
		for (int i = 0; i < gpsnum; i++) {
			GPSS[i].RO = sqrt((GPSS[i].xyz.x - Station.x) * (GPSS[i].xyz.x - Station.x) + (GPSS[i].xyz.y - Station.y) * (GPSS[i].xyz.y - Station.y) + (GPSS[i].xyz.z - Station.z) * (GPSS[i].xyz.z - Station.z));
		}
		for (int i = gpsnum; i < nums; i++) {
			BDSS[i - gpsnum].RO = sqrt((BDSS[i - gpsnum].xyz.x - Station.x) * (BDSS[i - gpsnum].xyz.x - Station.x) + (BDSS[i - gpsnum].xyz.y - Station.y) * (BDSS[i - gpsnum].xyz.y - Station.y) + (BDSS[i - gpsnum].xyz.z - Station.z) * (BDSS[i - gpsnum].xyz.z - Station.z));
		}
		//设置OMC向量矩阵W S*1
		double w;
		for (int i = 0; i < gpsnum; i++) {
			w = ((GPSS[i].xyz.x - Station.x) * (GPSS[i].xyzDot.x - StationSpeed.x) + (GPSS[i].xyz.y - Station.y) * (GPSS[i].xyzDot.y - StationSpeed.y) + (GPSS[i].xyz.z - Station.z) * (GPSS[i].xyzDot.z - StationSpeed.z)) / GPSS[i].RO + l_[i] * GPSS[i].xyzDot.x + m_[i] * GPSS[i].xyzDot.y + n_[i] * GPSS[i].xyzDot.z + c * GPSS[i].DEL_cld;
			W_.push_back(w);
		}
		for (int i = gpsnum; i < nums; i++) {
			w= ((BDSS[i - gpsnum].xyz.x - Station.x) * (BDSS[i - gpsnum].xyzDot.x - StationSpeed.x) + (BDSS[i - gpsnum].xyz.y - Station.y) * (BDSS[i - gpsnum].xyzDot.y - StationSpeed.y) + (BDSS[i - gpsnum].xyz.z - Station.z) * (BDSS[i - gpsnum].xyzDot.z - StationSpeed.z)) / BDSS[i - gpsnum].RO + l_[i - gpsnum] * BDSS[i - gpsnum].xyzDot.x + m_[i - gpsnum] * BDSS[i - gpsnum].xyzDot.y + n_[i - gpsnum] * BDSS[i - gpsnum].xyzDot.z + c * BDSS[i - gpsnum].DEL_cld;
			W_.push_back(w);
		}
		W.setMatrix(nums, 1, W_.data());
		//最小二乘
		X = reweightedLeastSquares(B, W, P);
		//接收机速度,钟漂数值更新,无须迭代计算
		StationSpeed.x += DataNormalize(X.data[0][0]);
		StationSpeed.y += DataNormalize(X.data[1][0]);
		StationSpeed.z += DataNormalize(X.data[2][0]);
		bd_gps += DataNormalize(X.data[3][0]);
		bd_bds += DataNormalize(X.data[4][0]);
		//精度评定
		GDOP = sqrt(Q.trace());
		PDOP = sqrt(Q.trace() - Q.data[3][3] * Q.data[3][3]);
		//输出接收机速度、钟漂及其精度信息
		cout << "接收机速度:" << StationSpeed.x << "," << StationSpeed.y << "," << StationSpeed.z << endl;
		cout << "接收机钟漂 " << "GPS系统:" << bd_gps << " 北斗系统:" << bd_bds << endl;
		cout << "精度信息 GDOP:" << GDOP << " PDOP:" << PDOP << endl;
	}

	void SinglePointPosition(string s) {
		SelectSatellite(s);
		SatelliteCount(sel);
		if (isSinglePosition(sel)) {
			cout << "单系统单点定位：" << endl;
			DATA_INIT_SSP();
			GetSSPP();
			GetSSPST();
		}
		else {
			cout << "组合系统单点定位：" << endl;
			DATA_INIT_CSP();
			GetCSPP();
			GetCSPST();
		}
	}
	//文件读取
	void FileRead_GMN(const char* FileName) {
		ifstream file(FileName);
		string line, s1, s2, s3, temp;
		int lineNumber = 0, cot = -1;
		int hdr_flag = 1;
		int flag = 0;
		int i = 0;
		if (!file.fail()) {
			while (getline(file, line)) {
				lineNumber++;
				if (cot > 2880)break;
				s1 = line.substr(60);
				s3 = line.substr(60, 20);
				if (hdr_flag == 1) {
					//信息头
					if (s1 == "IONOSPHERIC CORR    ") {
						s2 = line.substr(0, 4);
						for (int i_ = 0; i_ < 4; i_++) {
							if (line.substr(6 + 12 * i_, 11) == "           ") {
								corr.TEMP.push_back(0);
							}
							else
								corr.TEMP.push_back(my_stod(line.substr(6 + 12 * i_, 11)));

						}
						if (corr.TEMP.size() >= 4) {
							if (s2 == "GPSA") { corr.GPSAlpha_n[0] = corr.TEMP[0], corr.GPSAlpha_n[1] = corr.TEMP[1], corr.GPSAlpha_n[2] = corr.TEMP[2], corr.GPSAlpha_n[3] = corr.TEMP[3]; }
							else if (s2 == "GPSB") { corr.GPSBeta_n[0] = corr.TEMP[0], corr.GPSBeta_n[1] = corr.TEMP[1], corr.GPSBeta_n[2] = corr.TEMP[2], corr.GPSBeta_n[3] = corr.TEMP[3]; }
							else if (s2 == "BDSA") { corr.BDSAlpha_n[0] = corr.TEMP[0], corr.BDSAlpha_n[1] = corr.TEMP[1], corr.BDSAlpha_n[2] = corr.TEMP[2], corr.BDSAlpha_n[3] = corr.TEMP[3]; }
							else if (s2 == "BDSB") { corr.BDSBeta_n[0] = corr.TEMP[0], corr.BDSBeta_n[1] = corr.TEMP[1], corr.BDSBeta_n[2] = corr.TEMP[2], corr.BDSBeta_n[3] = corr.TEMP[3]; }
							corr.TEMP.clear();
						}
						else {
							cout << "corr size is less than 4!" << endl;
						}
						continue;
					}

				}
				if (s3 == "END OF HEADER       ") {
					hdr_flag = 0;
					s3.clear();
					continue;
				}
				if (hdr_flag == 0) {
					if ((line[0]=='G'|| line[0] == 'C') && isdigit(line[1]) && isdigit(line[2])) {
						flag = 1;
						if (i == 0&&flag==1) {
							cot++;
							sat[cot].gmn.PRN = line.substr(0, 3);
							sat[cot].gmn.TOC.year = stoi(line.substr(4, 4));
							sat[cot].gmn.TOC.month = stoi(line.substr(9, 2));
							sat[cot].gmn.TOC.day = stoi(line.substr(12, 2));
							sat[cot].gmn.TOC.hour = stoi(line.substr(15, 2));
							sat[cot].gmn.TOC.minute = stoi(line.substr(18, 2));
							sat[cot].gmn.TOC.second = my_stod(line.substr(21, 2));
							sat[cot].gmn.ClkBias = my_stod(line.substr(23, 19));
							sat[cot].gmn.ClkDrift = my_stod(line.substr(42, 19));
							sat[cot].gmn.ClkDriftRate = my_stod(line.substr(61, 19));
							i++;
							continue;
						}
					}
					
						 if (i == 1 && flag == 1) {
							sat[cot].gmn.IODE = my_stod(line.substr(4, 19));
							sat[cot].gmn.Crs = my_stod(line.substr(23, 19));
							sat[cot].gmn.DeltaN = my_stod(line.substr(42, 19));
							sat[cot].gmn.M0 = my_stod(line.substr(61, 19));
							i++;
							continue;
						}
						 if (i == 2 && flag == 1) {
							sat[cot].gmn.Cuc = my_stod(line.substr(4, 19));
							sat[cot].gmn.e = my_stod(line.substr(23, 19));
							sat[cot].gmn.Cus = my_stod(line.substr(42, 19));
							sat[cot].gmn.SqrtA = my_stod(line.substr(61, 19));
							i++;
							continue;
						}
						 if (i == 3 && flag == 1) {
							sat[cot].gmn.TOE = my_stod(line.substr(4, 19));
							sat[cot].gmn.Cic = my_stod(line.substr(23, 19));
							sat[cot].gmn.Omega = my_stod(line.substr(42, 19));
							sat[cot].gmn.Cis = my_stod(line.substr(61, 19));
							i++;
							continue;
						}
						 if (i == 4 && flag == 1) {
							sat[cot].gmn.i0 = my_stod(line.substr(4, 19));
							sat[cot].gmn.Crc = my_stod(line.substr(23, 19));
							sat[cot].gmn.omega = my_stod(line.substr(42, 19));
							sat[cot].gmn.OmegaDot = my_stod(line.substr(61, 19));
							i++;
							continue;
						}
						 if (i == 5 && flag == 1) {
							sat[cot].gmn.iDot = my_stod(line.substr(4, 19));
							sat[cot].gmn.CodesOnL2Channel = my_stod(line.substr(23, 19));
							sat[cot].gmn.GPSWeek = my_stod(line.substr(42, 19));
							sat[cot].gmn.L2PDataFlag = my_stod(line.substr(61, 19));
							i++;
							continue;
						}
						 if (i == 6 && flag == 1) {
							sat[cot].gmn.SVAccuracy = my_stod(line.substr(4, 19));
							sat[cot].gmn.SVHealth = my_stod(line.substr(23, 19));
							sat[cot].gmn.TGD = my_stod(line.substr(42, 19));
							sat[cot].gmn.IODC = my_stod(line.substr(61, 19));
							i++;
							continue;
						}
						 if (i == 7 && flag == 1) {
							sat[cot].gmn.TransTimeOfMsg = my_stod(line.substr(4, 19));
							sat[cot].gmn.Spare1 = my_stod(line.substr(23, 19));
							i = 0;
							flag = 0;
							continue;
						}
				}
				
			}
			file.close();
			cout << "文件已读取,请选择卫星:" << endl;
		}
		else { cout << "无法打开文件!" << endl; }
		
	}
};

