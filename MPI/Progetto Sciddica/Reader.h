#ifndef READER_H
#define READER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

class Reader
{
public:

    Reader(const char * path) :path(path), data(NULL)
    {


    }

    ~Reader()
    {
        if (data != NULL)
            delete [] data;
    }


    void loadFromFile ()
    {

        string line;
        ifstream myfile (path);
        int i = 0;

        int index=0;

        if (myfile.is_open())
        {
            while ( getline (myfile,line) )
            {

                vector<string> sep = split(line, ' ');
                switch (i) {
                case 0:
                    nCols= std::stoi(sep[1]);
                    break;
                case 1:
                    nRows= std::stoi(sep[1]);
                    data = new float [nRows*nCols];
                    break;
                default:

                    for(int col=0; col<nCols; col++)
                    {
                        data[index] = std::stof(sep[col]);
                        index++;
                    }

                    break;
                }

                i++;
            }
            myfile.close();
        }

    }

    float* getData ()
    {
        return data;
    }


    friend std::ostream & operator <<( std::ostream &os, const Reader &reader )
    {

        os<<"Matrix di size: "<<reader.nRows<<" X "<<reader.nCols<<"\n";

        for (int i = 0; i < reader.nRows; ++i) {
            for (int j = 0; j < reader.nCols; ++j) {
                os<<reader.data[i*reader.nCols+ j]<<" ";


            }
            os<<"\n";
        }
        return os;
    }




    unsigned int nRows;
    unsigned int nCols;
protected:
    float * data;
    const char* path;


    // You could also take an existing vector as a parameter.
    vector<string> split(string str, char delimiter) {
        vector<string> internal;
        stringstream ss(str); // Turn the string into a stream.
        string tok;

        while(getline(ss, tok, delimiter)) {

            internal.push_back(tok);
        }

        return internal;
    }

};

#endif
