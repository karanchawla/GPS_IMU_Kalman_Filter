//
//  paramReader.hpp
//  Fusion
//
//  Created by Karan on 6/6/18.
//  Copyright � 2018 Karan. All rights reserved.

#ifndef paramReader_hpp
#define paramReader_hpp

#endif /* paramReader_hpp */

#include <fstream>
#include <map>
#include <vector>

using namespace std;

class ParameterReader
{
public:
    ParameterReader( string filename="parameters.txt" )
    {
        ifstream fin( filename.c_str() );
        if (!fin)
        {
            cerr<<"parameter file does not exist."<<endl;
            return;
        }
        while(!fin.eof())
        {
		   cout << "------ Reading in Parameter File...\r\n";
            string str;
            getline( fin, str );
		  cout << "	Line Read: " << str << endl;
            if (str[0] == '#')
            {
                continue;
            }

            int pos = str.find("=");
            if (pos == -1){
		  	cout << "pos found = -1 ---- Continuing loop...\r\n";
                continue;
		 }
            string key = str.substr( 0, pos );
            string value = str.substr( pos+1, str.length() );

            this->data[key] = value;

		  cout << "	Key Found with Value: " << key << " -> " << value << endl;
		  cout << "	Stored data mapping:   key (" << key << ") ------- value(" << this->data[key] << ")\r\n";

            if ( !fin.good() ){
			  cout<<"\r\n";
			  break;
		  }

        }
    }
    string getData( string key )
    {
        map<string, string>::iterator iter;
	   iter = this->data.find(key.c_str());
	   std::cout << "Searching for key (" << key.c_str() << ") => " << this->data[key] << '\n';
        if (iter == this->data.end())
        {
            cerr<<"	Parameter name "<< key <<" not found!"<<endl;
            return string("NOT_FOUND");
        }
        return iter->second;
    }
public:
    map<string, string> data;
};
