
#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "Matrix.h"
// GLM Mathematics
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <climits>

class Temperature : public Matrix
{

public:

    Temperature (const char * path) :Matrix (path)
    {


        cout<<"max = "<<maximum << "min= "<<minimum<<endl;
        setMin();
        setMax();


        size = (coords.nRows+1) * (coords.nCols+1);
        temperatureColor = new float [size];

        computeColor();

    }

    void printColor ()
    {
        int globalIndex = 0;
        for (int i = 0; i < coords.nRows; ++i) {
            for (int j = 0; j < coords.nCols; ++j) {
                cout<<temperatureColor[globalIndex++]<<" ";
                if (globalIndex%coords.nCols ==0)
                    cout<<endl;

            }
        }

    }

    float* getColors()
    {
        return temperatureColor;
    }

    ~Temperature()
    {
        delete [] temperatureColor;
    }


    size_t getSize ()
    {
        return size;
    }



    void updateColor (float* temperature)
    {
        //x' = x - min / (max-min)


//        for (int i = 0; i < size; ++i) {
//            if(temperatureColor[i] == 0.0f)
//                temperatureColor[i]=0.5f;
//            else
//                temperatureColor[i] = 0.0f;
//        }
        updateMinMax(temperature);
        int globalIndex = 0;
        for (int i = 0; i < coords.nRows+1; ++i) {
            for (int j = 0; j < coords.nCols+1; ++j) {
                if (i>=coords.nRows || j>= coords.nCols || temperature[globalIndex] == 0.0f)
                {
                    temperatureColor[globalIndex++] = 0.0f;
                }
                else
                    temperatureColor[globalIndex++] = (temperature[globalIndex] - minimum) / (maximum - minimum);
            }
        }
    }




protected:
    float* temperatureColor;

    float minimum;
    float maximum;

    size_t size;

    void setMin ()
    {
        bool first = false;
        for (int i = 0; i < coords.nRows; ++i) {
            for (int j = 0; j < coords.nCols; ++j) {
                if (data[i][j] != NODATA_value)
                {
                    if (!first)
                    {
                        minimum = data[i][j];
                        first = true;
                    }
                    if (data[i][j] < minimum)
                    {
                        minimum = data[i][j];
                    }
                }
            }
        }
    }

    void setMax ()
    {
        bool first = false;
        for (int i = 0; i < coords.nRows; ++i) {
            for (int j = 0; j < coords.nCols; ++j) {
                if (data[i][j] != NODATA_value)
                {
                    if (!first)
                    {
                        maximum = data[i][j];
                        first = true;
                    }
                    if (data[i][j] > maximum)
                    {
                        maximum = data[i][j];
                    }
                }
            }
        }
    }

    void computeColor()
    {
        //x' = x - min / (max-min)
        int globalIndex = 0;
        for (int i = 0; i < coords.nRows+1; ++i) {
            for (int j = 0; j < coords.nCols+1; ++j) {
                if (i>=coords.nRows || j>= coords.nCols || data[i][j] == 0.0f)
                {
                    temperatureColor[globalIndex++] = 0.0f;
                }
                else
                    temperatureColor[globalIndex++] = (data[i][j] - minimum) / (maximum - minimum);
            }
        }

    }


    void updateMinMax (float* temperature)
    {

        bool first = false;
        for (int i = 0; i < coords.nRows*coords.nCols; ++i) {

                if (temperature[i] != NODATA_value)
                {
                    if (!first)
                    {
                        minimum = temperature[i];
                        maximum = temperature[i];
                        first = true;
                    }
                    if (temperature[i] < minimum)
                    {
                        minimum = temperature[i];
                    }

                    if (temperature[i] > maximum)
                    {
                        maximum = temperature[i];
                    }
                }
            }

    }






};

#endif
