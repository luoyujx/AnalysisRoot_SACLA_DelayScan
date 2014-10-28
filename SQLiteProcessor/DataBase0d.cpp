#include "DataBase0d.h"
#include <sstream>
//#include <stdlib.h>
#include <algorithm>

DataBase0d::DataBase0d(void)
{
}
DataBase0d::~DataBase0d(void)
{
}
//Open SQLite data base
void DataBase0d::Open(const string &dbFileName)
{
	int result = sqlite3_open_v2(dbFileName.c_str(), &dataBase, SQLITE_OPEN_READONLY, 0 );
	if (result == SQLITE_OK) 
	{
		std::cout << dbFileName << " has been opened." << std::endl;
	}
	else
	{
		std::cout << "Error!!: " <<::sqlite3_errmsg(dataBase) << std::endl;
	}
}
//
void DataBase0d::LoadData(unsigned int firstTag, unsigned int lastTag, vector<string> &fields)
{
	//copy field names to member
	fieldNames.resize(fields.size());
	std::copy(fields.begin(), fields.end(), fieldNames.begin());
	//making SQL statement
	std::stringstream command;
	command << "select tag";
	for (unsigned int i = 0; i < fields.size(); i++)
	{
		command << ", [" << fields[i] << "] ";
	}
	command << "from bldata where tag between ";
	command << firstTag << " and " << lastTag << ";";
	std::cout << "SQL command: " << command.str() << std::endl;
	//compile SQL command
	sqlite3_prepare_v2(dataBase, command.str().c_str(), -1, &statement, 0);
	//
	int tag = 0;
	std::vector<double> fieldData(fields.size());
	while (sqlite3_step(statement) == SQLITE_ROW) 
	{
		tag = sqlite3_column_int(statement,0);
		for (unsigned int i = 0; i < fields.size(); i++)
		{
			const char* dataTmp = reinterpret_cast<const char*>(sqlite3_column_text(statement, i+1));
			char *pEnd = NULL;
			fieldData[i] = strtod(dataTmp, &pEnd);
			if (!pEnd)
			{
				fieldData[i] = std::numeric_limits<double>::quiet_NaN();
			}
		}		
		//Load data to internal table
		table.insert(std::pair<unsigned int, std::vector<double> >(tag, fieldData));
	}
	//Close SQLite Data base
	sqlite3_finalize(statement);
	sqlite3_close(dataBase);
}
//
void DataBase0d::ShowTable()
{
	std::cout << "***Show table***" <<std::endl;
	std::cout << "Tag\t";
	std::for_each(fieldNames.begin(), fieldNames.end(), outputTabDelim(std::cout, fieldNames.size()));
	for (tableIntDouble::iterator it=table.begin(); it != table.end(); it++)
	{
		std::cout << it->first << "\t";
		std::for_each(it->second.begin(), it->second.end(),outputTabDelim(std::cout, it->second.size()));
	}
}
//
double DataBase0d::GetData(unsigned int tag, unsigned int fieldNumber)
{
	tableIntDouble::iterator itFind = table.find(tag);
	if (itFind == table.end())
	{
		return std::numeric_limits<double>::quiet_NaN();
	}
	return itFind->second[fieldNumber];
}
