
#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "../Reader.h"
// GLM Mathematics
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <climits>



class Temperature : public Reader
{

public:

    Temperature (const char * path) :Reader (path)
    {

        loadFromFile();
        fillMatrix();

        setMin();
        setMax();

        cout<<"max = "<<maximum << "min= "<<minimum<<endl;

        coords.nCols = nCols+1;
        coords.nRows = nRows+1;

        size = coords.nCols * coords.nRows;
        temperatureColor = new double [size];

//        updateColor(this->getDataLinear());
        computeColor();



    }

    void printColor ()
    {
        int globalIndex = 0;
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                cout<<temperatureColor[globalIndex++]<<" ";
//                if (globalIndex%nCols ==0)

            }
            cout<<endl;
        }

    }

    double* getColors()
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



    void updateColor (double* temperature)
    {
        //x' = x - min / (max-min)


        updateMinMax(temperature);
        int globalIndex = 0;
        for (int i = 0; i < coords.nRows; ++i) {
            for (int j = 0; j < coords.nCols; ++j) {
                if (i>=nRows || j>= nCols || temperature[i*nCols+j] == 0.0f)
                {
                    temperatureColor[i*coords.nCols+j] = 0.0f;
                }
                else
                {

                    temperatureColor[i*coords.nCols+j] = (temperature[i*nCols+j] - minimum) / (maximum - minimum);
                }
            }
        }
    }




protected:
    double* temperatureColor;

    double minimum;
    double maximum;

    Coordinates coords;
    size_t size;

    void setMin ()
    {
        bool first = false;
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
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
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
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

        double * tmp = this->getDataLinear();
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                if (tmp[i*nCols+j] != data[i][j])
                    cout<<"Mannaia la puttana "<<i<<"  "<<j<<endl;
            }
        }

//        double* tmp = getDataLinear();
        for (int i = 0; i < coords.nRows; ++i) {
            for (int j = 0; j < coords.nCols; ++j) {
                if (i>=nRows || j>= nCols || data[i][j] == 0.0f)
                {
                    temperatureColor[i*coords.nCols+j] = 0.0f;
                }
                else
                    temperatureColor[i*coords.nCols+j] = (data[i][j] - minimum) / (maximum - minimum);
            }
        }

    }


    void updateMinMax (double* temperature)
    {

        bool first = false;
        for (int i = 0; i < nRows*nCols; ++i) {

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
