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
	//mysql_close(dataBase);
}
//Open SQLite data base
void DataBase0d::Open(const string &dbFileName)
{
	int result = sqlite3_open_v2(dbFileName.c_str(), &dataBaseL, SQLITE_OPEN_READONLY, 0 );
	if (result == SQLITE_OK) 
	{
		std::cout << dbFileName << " has been opened." << std::endl;
	}
	else
	{
		std::cout << "Error!!: " <<::sqlite3_errmsg(dataBaseL) << std::endl;
	}
}

void DataBase0d::Open(const char *dbFileName)
{
	int result = sqlite3_open_v2(dbFileName, &dataBaseL, SQLITE_OPEN_READONLY, 0 );
	if (result == SQLITE_OK) 
	{
		std::cout << dbFileName << " has been opened." << std::endl;
	}
	else
	{
		std::cout << "Error!!: " <<::sqlite3_errmsg(dataBaseL) << std::endl;
	}
}
void DataBase0d::Connect(const char* DBHOST, const char* DBUSER, const char* DBPASS, const char* DBNAME)
{	
	dataBaseM = mysql_init(NULL);
	if (!mysql_real_connect(dataBaseM, DBHOST, DBUSER, DBPASS, DBNAME, 3306, NULL, 0))
	{
		std::cout << "Error!!: " << mysql_error(dataBaseM) << std::endl;
	}
	else
	{
		std::cout << "MySQL server: " << DBHOST << " Connection succeeded." << std::endl;
	}
}
//
void DataBase0d::LoadDataL(unsigned int firstTag, unsigned int lastTag, vector<string> &fields)
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
	sqlite3_prepare_v2(dataBaseL, command.str().c_str(), -1, &statement, 0);
	//
	int tag = 0;
	std::vector<double> fieldData(fields.size());
	while (sqlite3_step(statement) == SQLITE_ROW) 
	{
		tag = sqlite3_column_int(statement,0);
		for (unsigned int i = 0; i < fields.size(); i++)
		{
			const char* dataTmp = reinterpret_cast<const char*>(sqlite3_column_text(statement, i+1));
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
	//Close SQLite Data base
	sqlite3_finalize(statement);
	sqlite3_close(dataBaseL);
}
void DataBase0d::LoadDataM(unsigned int runNum, vector<string> &fields)
{
	//making SQL statement
	std::stringstream command;
	command << "select StartTag, EndTag ";
	command << "from RunInfo where run = ";
	command << runNum << ";";
	std::cout << "SQL command: " << command.str() << std::endl;
	//Query SQL command
	int result = mysql_query(dataBaseM, command.str().c_str());
	//Get result
	res = mysql_use_result(dataBaseM);
	if ((row = mysql_fetch_row(res)) != NULL)
	{
		unsigned int firstTag = atoi(row[0]);
		unsigned int lastTag = atoi(row[1]);
		mysql_free_result(res);
		LoadDataM(firstTag, lastTag, fields);
	}
	
}
void DataBase0d::LoadDataM(unsigned int firstTag, unsigned int lastTag, vector<string> &fields)
{
	//copy field names to member
	fieldNames.resize(fields.size());
	std::copy(fields.begin(), fields.end(), fieldNames.begin());
	//making SQL statement
	std::stringstream command;
	command << "select tag";
	for (unsigned int i = 0; i < fields.size(); i++)
	{
		command << ", `" << fields[i] << "` ";
	}
	//command << "from bldata where tag between ";
	command << "from bldata where tag between ";
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
	mysql_close(dataBaseM);
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
//
std::pair<int,double> DataBase0d::GetStatusAndData(unsigned int tag, unsigned int fieldNumber)
{
	std::pair<int,double> retval(0,std::numeric_limits<double>::quiet_NaN());
	tableIntDouble::iterator itFind = table.find(tag);
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
