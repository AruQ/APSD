#ifndef READER_H
#define READER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;
struct Coordinates
{
    unsigned int nRows;
    unsigned int nCols;

};
class Reader
{
public:

    unsigned int nRows;
    unsigned int nCols;

    Reader(const char * path) :path(path), data(NULL)
    {


    }

    ~Reader()
    {
        for (int i = 0; i < nRows; ++i) {
            delete [] data[i];
        }
        delete [] data;
    }


    void loadFromFile ()
    {

        string line;
        ifstream myfile (path);
        int i = 0;
        int row=0;

        if (myfile.is_open())
        {
            while ( getline (myfile,line) )
            {

                vector<string> sep = split(line, ' ');
                switch (i) {
                case 0:
                    nCols= atoi(sep[1].c_str());
                    break;
                case 1:
                    nRows= atof(sep[1].c_str());
                    data = new float* [nRows];
                    break;
                case 2:
                    xllcorner= atof(sep[1].c_str());
                    break;
                case 3:
                    yllcorner= atof(sep[1].c_str());
                    break;
                case 4:
                    cellsize=atof(sep[1].c_str());
                    break;
                case 5:
                    NODATA_value= atof(sep[1].c_str());
                    break;
                default:
                    data[row] = new float[nCols];

                    for(int col=0; col<nCols; col++)
                    {
                        data[row][col] = atof(sep[col].c_str());
                    }
                    row++;
                    break;
                }

                i++;
            }
            myfile.close();
        }

    }



    float* getDataLinear ()
    {
        int globalIndex = 0;
        float * _data = new float[nRows* nCols];
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                _data[globalIndex++] = this->data[i][j];
            }
        }
        return _data;
    }


    void fillMatrix ()
    {
        if (nCols == nRows)
            return;
        float** tmp;
        if (nCols > nRows)
        {
            tmp = new float* [nCols];
            for (int i = 0; i < nCols; ++i) {
               tmp[i] = new float[nCols];
            }

            int diff = nCols-nRows;

            for(int i = 0;i <diff/2;i++)
                for(int j = 0;j< nCols;j++)
                    tmp[i][j] = data[0][j];

            for(int i = nRows+diff/2;i <nCols;i++)
                for(int j = 0;j< nCols;j++)
                    tmp[i][j] = data[nRows-1][j];


            for(int i = diff/2;i<nRows+diff/2;i++)
                for(int j = 0;j< nCols;j++)
                    tmp[i][j] = data[i-diff/2][j];


        }

        else
        {
            tmp = new float* [nRows];
            for (int i = 0; i < nRows; ++i) {
               tmp[i] = new float[nRows];
            }

            int diff = nRows-nCols;

            for(int i = 0;i <nRows;i++)
                for(int j = 0;j< diff/2;j++)
                    tmp[i][j] = data[i][0];

            for(int i = 0;i <nRows;i++)
                for(int j = nCols+diff/2;j< nRows;j++)
                    tmp[i][j] = data[i][nCols-1];


            for(int i = 0;i<nRows;i++)
                for(int j = diff/2;j< nCols+diff/2;j++)
                    tmp[i][j] = data[i][j-diff/2];

        }

        for(int i =0;i<nRows;i++)
        {
            delete [] data[i];
        }
        delete [] data;
        data = tmp;

        nCols = nCols > nRows? nCols:nRows;
        nRows = nCols > nRows? nCols:nRows;
    }

    void fillMatrixNODATAValue ()
    {
        if (nCols == nRows)
            return;
        float** tmp;
        if (nCols > nRows)
        {
            tmp = new float* [nCols];
            for (int i = 0; i < nCols; ++i) {
               tmp[i] = new float[nCols];
            }

            int diff = nCols-nRows;

            for(int i = 0;i <diff/2;i++)
                for(int j = 0;j< nCols;j++)
                    tmp[i][j] = NODATA_value;

            for(int i = nRows+diff/2;i <nCols;i++)
                for(int j = 0;j< nCols;j++)
                    tmp[i][j] = NODATA_value;


            for(int i = diff/2;i<nRows+diff/2;i++)
                for(int j = 0;j< nCols;j++)
                    tmp[i][j] = data[i-diff/2][j];


        }

        else
        {
            tmp = new float* [nRows];
            for (int i = 0; i < nRows; ++i) {
               tmp[i] = new float[nRows];
            }

            int diff = nRows-nCols;

            for(int i = 0;i <nRows;i++)
                for(int j = 0;j< diff/2;j++)
                    tmp[i][j] =NODATA_value;

            for(int i = 0;i <nRows;i++)
                for(int j = nCols+diff/2;j< nRows;j++)
                    tmp[i][j] = NODATA_value;


            for(int i = 0;i<nRows;i++)
                for(int j = diff/2;j< nCols+diff/2;j++)
                    tmp[i][j] = data[i][j-diff/2];

        }

        for(int i =0;i<nRows;i++)
        {
            delete [] data[i];
        }
        delete [] data;
        data = tmp;

        nCols = nCols > nRows? nCols:nRows;
        nRows = nCols > nRows? nCols:nRows;
    }

    float ** getData()
    {
        return data;
    }

    float getCellSize()
    {
        return cellsize;
    }

    friend std::ostream & operator <<( std::ostream &os, const Reader &matrix )
    {

        os<<"Matrix di size: "<<matrix.nRows<<" X "<<matrix.nCols<<"\n";

        for (int i = 0; i < matrix.nRows; ++i) {
            for (int j = 0; j < matrix.nCols; ++j) {
                os<<matrix.data[i][j]<<" ";


            }
            os<<"\n";
        }
        return os;
    }


protected:


    float xllcorner;
    float yllcorner;
    float cellsize;
    float NODATA_value;

    float ** data;
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
