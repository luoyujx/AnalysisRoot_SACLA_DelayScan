#include "DataBase0d.h"
#include "stdafx.h"
#include <sstream>
//#include <stdlib.h>
#include <algorithm>

DataBase0d::DataBase0d(void)
{

}
DataBase0d::~DataBase0d(void)
{
	//mysql_close(dataBaseM);
}
// Connect MySQL data base
void DataBase0d::Open(string fileName)
{
	dataFs.open(fileName, std::ios::in);
	if (dataFs.fail())
	{
		std::cout << "Can not open db text file." << fileName << std::endl;
	}
}
void DataBase0d::LoadData(int numOfFields)
{
	unsigned int tag = 0;
	std::vector<double> fieldData(numOfFields);
	char tmp[256];
	while (!dataFs.eof())
	{
		dataFs >> tag;
		//std::cout << tag << std::endl;
		for (int i = 0; i < numOfFields; i++)
		{
			dataFs >> fieldData[i];
		}
		dataFs.getline(tmp, 256);
		//Load data to internal table
		table.insert(std::pair<unsigned int, std::vector<double> >(tag, fieldData));
	}
}
//
void DataBase0d::ShowTable()
{
	std::cout << "***Show table***" << std::endl;
	//std::cout << "Tag\t";
	//std::for_each(fieldNames.begin(), fieldNames.end(), outputTabDelim(std::cout, fieldNames.size()));
	for (tableIntDouble::iterator it = table.begin(); it != table.end(); it++)
	{
		std::cout << it->first << "\t";
		std::for_each(it->second.begin(), it->second.end(), outputTabDelim(std::cout, it->second.size()));
	}
}
//
double DataBase0d::GetData(unsigned int tag, unsigned int fieldNumber)
{
	tableIntDouble::iterator itFind = table.find(tag);
	//tableIntDouble::iterator itFind = table.upper_bound(tag);
	//itFind--;
	if (itFind == table.end())
	{
		return std::numeric_limits<double>::quiet_NaN();
	}
	return itFind->second[fieldNumber];
}
//
std::pair<int, double> DataBase0d::GetStatusAndData(unsigned int tag, unsigned int fieldNumber)
{
	std::pair<int, double> retval(0, std::numeric_limits<double>::quiet_NaN());
	tableIntDouble::iterator itFind = table.find(tag);
	//tableIntDouble::iterator itFind = table.upper_bound(tag);
	//itFind--;
	if (itFind == table.end())
	{
		retval.first = 0;
		retval.second = std::numeric_limits<double>::quiet_NaN();
		return retval;
	}
	retval.first = 1;
	retval.second = itFind->second[fieldNumber];
	return retval;
}
//
//int DataBase0d::GetLatestTag(const char* TABLE)
//{
//	std::stringstream command;
//	command << "select max(Tag) from " << TABLE;
//	int result = mysql_query(dataBaseM, command.str().c_str());
//	res = mysql_use_result(dataBaseM);
//	int ltag;
//	while ((row = mysql_fetch_row(res)) != NULL)
//	{
//		ltag = atoi(row[0]);
//	}
//	mysql_free_result(res);
//	//mysql_close(dataBaseM);
//	return ltag - 1000;
//}
//
void DataBase0d::Close()
{
	dataFs.close();
}