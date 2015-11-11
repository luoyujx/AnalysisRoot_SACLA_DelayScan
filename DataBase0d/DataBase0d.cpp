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
void DataBase0d::Connect(const char* DBHOST, const char* DBUSER, const char* DBPASS, const char* DBNAME)
{
	dataBaseM = mysql_init(NULL);
	if (!mysql_real_connect(dataBaseM, DBHOST, DBUSER, DBPASS, DBNAME, 3306, NULL, 0))
	{
		std::cout << "Error!!: " << mysql_error(dataBaseM) << std::endl;
	}
	else
	{
		//std::cout << "MySQL server: " << DBHOST << " Connection succeeded." << std::endl;
	}
}
void DataBase0d::LoadDataM(int firstTag, int lastTag, vector<string> &fields, const char* TABLE)
{
	//copy field names to member
	fieldNames.resize(fields.size());
	std::copy(fields.begin(), fields.end(), fieldNames.begin());
	//making SQL statement
	std::stringstream command;
	command << "select Tag";
	for (unsigned int i = 0; i < fields.size(); i++)
	{
		command << ", `" << fields[i] << "` ";
	}
	command << "from " << TABLE << " where tag between ";
	command << firstTag << " and " << lastTag << ";";
	std::cout << "SQL command: " << command.str() << std::endl;
	//Query SQL command
	int result = mysql_query(dataBaseM, command.str().c_str());
	//Get result
	res = mysql_use_result(dataBaseM);
	int tag = 0;
	std::vector<double> fieldData(fields.size());
	while ((row = mysql_fetch_row(res)) != NULL)
	{
		tag = atoi(row[0]);
		//std::cout << tag << std::endl;
		for (unsigned int i = 0; i < fields.size(); i++)
		{
			const char* dataTmp = row[i + 1];
			if (dataTmp == NULL)
			{
				fieldData[i] = std::numeric_limits<double>::quiet_NaN();
			}
			else
			{
				char *pEnd = NULL;
				fieldData[i] = strtod(dataTmp, &pEnd);
				if (!pEnd)
				{
					fieldData[i] = std::numeric_limits<double>::quiet_NaN();
				}
			}
		}
		//Load data to internal table
		table.insert(std::pair<unsigned int, std::vector<double> >(tag, fieldData));
	}	
	//Close MySQL connection
	mysql_free_result(res);
	//mysql_close(dataBaseM);
}
//
void DataBase0d::ShowTable()
{
	std::cout << "***Show table***" << std::endl;
	std::cout << "Tag\t";
	std::for_each(fieldNames.begin(), fieldNames.end(), outputTabDelim(std::cout, fieldNames.size()));
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
int DataBase0d::GetLatestTag(const char* TABLE)
{
	std::stringstream command;
	command << "select max(Tag) from " << TABLE;
	int result = mysql_query(dataBaseM, command.str().c_str());
	res = mysql_use_result(dataBaseM);
	int ltag;
	while ((row = mysql_fetch_row(res)) != NULL)
	{
		ltag = atoi(row[0]);
	}
	mysql_free_result(res);
	//mysql_close(dataBaseM);
	return ltag - 1000;
}
//
void DataBase0d::CloseMySQL()
{
	mysql_close(dataBaseM);
}