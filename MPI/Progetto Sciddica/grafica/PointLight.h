#ifndef POINTLIGHT_H
#define POINTLIGHT_H



#include "Shader.h"


class Intensity {
private:
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;
public:
    Intensity ( glm::vec3 _ambient, glm::vec3 _diffuse, glm::vec3 _specular) :
        ambient(_ambient), diffuse(_diffuse), specular(_specular){}

    Intensity (const Intensity & _intensity)
    {
        this->ambient = _intensity.ambient;
        this->diffuse = _intensity.diffuse;
        this->specular = _intensity.specular;
    }

    void SetUniformData(Shader * shaderProgram, string lightVarName)
    {
        glUniform3f(glGetUniformLocation(shaderProgram->Program, (lightVarName+".ambient").c_str()),ambient.x, ambient.y, ambient.z);
        glUniform3f(glGetUniformLocation(shaderProgram->Program, (lightVarName+".diffuse").c_str()),diffuse.x, diffuse.y, diffuse.z);
        glUniform3f(glGetUniformLocation(shaderProgram->Program, (lightVarName+".specular").c_str()),specular.x, specular.y, specular.z);

    }


};

class PointLight
{
public:


    PointLight(int id, glm::vec4 _position, float _constant, float _linear, float _quadratic, Intensity _intensity):
        position(_position), constant(_constant), linear(_linear), quadratic (_quadratic), intensity(_intensity)

    {}

    void SetUniformData(Shader* shaderProgram, string lightVarName)
    {


//        std::cout<<shaderProgram->Program<<std::endl;
        glUniform3f(glGetUniformLocation(shaderProgram->Program, (lightVarName+".position").c_str()),position.x, position.y, position.z);

//        glUniform3f(glGetUniformLocation(shaderProgram->Program, (lightVarName+".viewPos").c_str()),position.x, position.y, position.z);
//        glUniform3f(glGetUniformLocation(shaderProgram->Program, (lightVarName+".color").c_str()),color.x, color.y, color.z);
        glUniform1f(glGetUniformLocation(shaderProgram->Program, (lightVarName+".costant").c_str()), constant);
        glUniform1f(glGetUniformLocation(shaderProgram->Program, (lightVarName+".linear").c_str()), linear);
        glUniform1f(glGetUniformLocation(shaderProgram->Program, (lightVarName+".quadratic").c_str()),quadratic);

        intensity.SetUniformData(shaderProgram, lightVarName);

    }

    void setPosition(glm::vec3 _position)
    {
        position.x = _position.x;
        position.y = _position.y;
        position.z = _position.z;
    }




    friend std::ostream & operator <<( std::ostream &os, const PointLight &pointLight )
      {
          os<<"PointLight position = "<<pointLight.position.x<< " "<<pointLight.position.y<<" "<<pointLight.position.z<<"\n";
          os<<"PointLight costant = "<<pointLight.constant<<"\n";
          os<<"PointLight linear = "<<pointLight.linear<<"\n";
          os<<"PointLight quadratic = "<<pointLight.quadratic<<"\n";
          return os;
      }

    glm::vec4& getPosition()
    {
        return position;
    }

private:

    glm::vec4 position;

    Intensity intensity;

    float constant;
    float linear;
    float quadratic;

    int id;

};


#endif
