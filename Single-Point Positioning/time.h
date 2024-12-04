#pragma once
#include<iostream>
#include<math.h>
#include<iomanip>
#include<string>

using namespace std;
/*儒略日定义 从公元前4713年1月1日正午开始的天数*/
class JulianDay {
public:

	unsigned long day;//具体的天数
	long sn;//当前天的秒数部分
	double tos;//当前时间的亚秒部分
	JulianDay(int a = 0, int b = 0, double c = 0) : day(a), sn(b), tos(c) {}

	double JulianDay_Output() {
		double julianday = this->day + (this->sn + this->tos) / 86400.0;
		return julianday;
	}
};
/*简化儒略日 从1858年11月17日子夜开始的天数*/
class MJulianDay {
public:

	unsigned long day;
	long sn;
	double tos;
	MJulianDay(int a = 0, int b = 0, double c = 0) : day(a), sn(b), tos(c) {}

	double MJulianDay_Output() {
		double mjulianday = this->day + (this->sn + this->tos) / 86400.0;
		return mjulianday;
	}
};

bool isLeapYear(unsigned short y) {
	return (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0);
}
int daysInMonth(int month, int year) {
	int days[] = { 31, isLeapYear(year) ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	return days[month - 1];
}

/* 通用时间定义*/
class  COMMONTIME{
public:
	unsigned short year;
	unsigned short month;
	unsigned short day;
	unsigned short hour;
	unsigned short minute;
	double second;
	COMMONTIME(int a = 0, int b = 0, int c = 0, int d = 0, int e = 0, double f = 0) :year(a), month(b), day(c), hour(d), minute(e), second(f) {}

	void COMTPrint() {
		std::cout << std::setw(4) << std::setfill('0') << year << "/"
			<< std::setw(2) << std::setfill('0') << month << "/"
			<< std::setw(2) << std::setfill('0') << day << " "
			<< std::setw(2) << std::setfill('0') << hour << ":"
			<< std::setw(2) << std::setfill('0') << minute << ":"
			<< std::fixed << std::setprecision(0) << std::setw(2) << (int)second << endl;
	}
	unsigned long totalDays() const {
		static const unsigned short daysPerMonth[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
		unsigned long days = 0;
		for (unsigned short m = 1; m < month; ++m) {
			days += daysPerMonth[m - 1];
			if (m == 2 && isLeapYear(year)) {
				days += 1; // 闰年二月加一天
			}
		}
		days += day - 1;

		return days;
	}
	double COM2GPSTime() {
		// GPS时间起点
		COMMONTIME gpsEpoch;
		gpsEpoch.year = 1980;
		gpsEpoch.month = 1;
		gpsEpoch.day = 6;
		unsigned long daysSinceGPSEpoch = this->totalDays() - gpsEpoch.totalDays();
		unsigned long gpsWeek = daysSinceGPSEpoch / 7;
		unsigned long daysInCurrentWeek = daysSinceGPSEpoch % 7;
		double gpsWeekSeconds = daysInCurrentWeek * 86400 + hour * 3600 + minute * 60 + second;
		return gpsWeekSeconds;
	}
	double calculateJulianDay() {
		int dayOfYear = day;
		for (int m = 1; m < month; ++m) dayOfYear += daysInMonth(m, year);
		return dayOfYear + (hour + minute / 60.0 + second / 3600.0) / 24.0;
	}
	double operator-(COMMONTIME& object) {
		unsigned long days1 = this->totalDays(); // 自1970年1月1日以来的天数
		unsigned long seconds1 = days1 * 86400 + this->hour * 3600 + this->minute * 60 + this->second;
		unsigned long days2 = object.totalDays();
		unsigned long seconds2 = days2 * 86400 + object.hour * 3600 + object.minute * 60 + object.second;
		return static_cast<double>(seconds1 - seconds2);
	}

};
COMMONTIME JulianDay2COM2024(double days) {
	COMMONTIME ct;
	ct.year = 2024;

	int D = static_cast<int>(std::floor(days)); // 整数部分即为天数
	double fractional = days - D; // 小数部分即为时间的亚秒部分

	// 计算具体日期
	static const int daysInMonth[12] = { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

	ct.month = 1; // 从1月开始计算
	while (D > daysInMonth[ct.month - 1]) {
		D -= daysInMonth[ct.month - 1]; // 扣除当前月的天数
		ct.month++; // 进入下一个月
	}
	ct.day = D; // 剩余天数即为当天的日期

	// 计算时分秒
	double totalSeconds = fractional * 86400; // 将小数部分转换为秒
	ct.hour = static_cast<unsigned short>(totalSeconds / 3600); // 计算小时
	ct.minute = static_cast<unsigned short>((totalSeconds - ct.hour * 3600) / 60); // 计算分钟
	ct.second = totalSeconds - ct.hour * 3600 - ct.minute * 60; // 计算秒

	return ct;
}

/*GPS时定义 起点是1980年1月6日00时00分00秒，用周和周内秒来表示*/
class GPSTime
{
public:
	int wn;
	long sn;
	double tos;
	GPSTime(unsigned short a = 0, int b = 0, double c = 0) :wn(a), sn(b), tos(c) {}
	
};
double frac(double x) {
	return x - int(x);
}
//儒略日转简化儒略日
double JulianDay2MJulianDay(double JD) {
	return JD - 2400000.5;
}
//简化儒略日转儒略日
double MJulianDay2JulianDay(double MJD) {
	return MJD + 2400000.5;
}
//1.通用时到儒略日的转换
void CommonTime2JulianDay(COMMONTIME ct, JulianDay& jd) {
	int y = 0, m = 0;
	y = (ct.month <= 2) ? (ct.year - 1) : ct.year;
	m = (ct.month <= 2) ? (ct.month + 12) : ct.month;
	double Day = int(365.25 * y) + int(30.6001 * (m + 1)) + ct.day + ct.hour / 24 + 1720981.5;//UT/24
	jd.day=int(Day);
	jd.sn = int(frac(Day) * 86400.0);
	jd.tos = frac(frac(Day) * 86400.0);
}
//2.儒略日转通用时
void JulianDay2CommonTime(COMMONTIME& ct, JulianDay& jd) {
	double x = jd.JulianDay_Output();
	int a = int(x + 0.5);
	int b = a + 1537;
	int c = int((b - 122.1) / 365.25);
	int d = int(365.25 * c);
	int e = int((b - d) / 30.6001);
	ct.day = b - d - int(30.6001 * e) + frac(x + 0.5);
	ct.month = e - 1 - 12 * int(e / 14);
	ct.year = c - 4715 - int((7 + ct.month) / 10);
	ct.hour = (jd.sn / 3600 + 12) % 24;
	ct.minute = (jd.sn % 3600) / 60;
	ct.second = jd.sn % 60 + jd.tos;
}
//3.GPS时转儒略日
void GPSTime2JulianDay(GPSTime gt, JulianDay& jd) {
	double MJD = 44244 + gt.wn * 7 + gt.sn / 86400;
	jd.day =int (MJulianDay2JulianDay(MJD));
	jd.sn = int(frac(MJulianDay2JulianDay(MJD)) * 86400.0);
	jd.tos = frac(frac(MJulianDay2JulianDay(MJD)) * 86400.0);
}
//4. 儒略日到GPS时
void JulianDay2GPSTime(JulianDay jd, GPSTime& gt) {
	double MJD = JulianDay2MJulianDay(jd.JulianDay_Output());
	gt.wn = int((MJD - 44244) / 7);
	gt.sn = (MJD - 44244 - gt.wn * 7) * 86400;
}
//5. 通用时到GPS时
void CommonTime2GPSTime(COMMONTIME ct, GPSTime& gt) {
	JulianDay jd;
	CommonTime2JulianDay(ct, jd);
	JulianDay2GPSTime(jd, gt);
}
//6. GPS时到通用时
void GPSTime2CommonTime(COMMONTIME& ct, GPSTime gt) {
	JulianDay jd;
	GPSTime2JulianDay(gt, jd);
	JulianDay2CommonTime(ct, jd);
}
//通用时转年积日