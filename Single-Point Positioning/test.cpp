#include"func.h"
int main() {
	DATA data;
	string s;
	data.FileRead_GMO("TWTF00TWN_R_20242090000_01D_30S_MO.rnx");
	data.FileRead_GMN("BRDM00DLR_S_20242090000_01D_MN.rnx");
	std::getline(std::cin, s);
	data.SinglePointPosition(s);
	return 0;
}
