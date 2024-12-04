#include"func.h"
int main() {
	DATA data;
	string s;
	data.FileRead_GMN("BRDM00DLR_S_20242090000_01D_MN.rnx");
	std::getline(std::cin, s);
	data.SinglePointPosition(s);
	system("pause");
	return 0;
}