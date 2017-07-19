#ifndef MYTEXTURE_H
#define MYTEXTURE_H

#include <GL/glew.h>

// Other Libs
#include <SOIL.h>
#include "Shader.h"

#include <string>

using namespace std;


class MyTexture
{
public:
    MyTexture(const char* _path)
    {
        id=instanceID;
        instanceID++;

        path =(char*)_path;
        glGenTextures(1, &texture);
        loadImage();
    }



    void setParameters (int _wrapSParameter, int _wrapTParameter, int _minFilter, int _magFilter)
    {
        wrapSParameter = _wrapSParameter;
        wrapTParameter = _wrapTParameter;
        minFilter = _minFilter;
        magFilter = _magFilter;

        glBindTexture(GL_TEXTURE_2D, texture);

        // Set our texture parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapSParameter);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapTParameter);
        // Set texture filtering
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,minFilter);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);



        glBindTexture(GL_TEXTURE_2D, 0);
    }



    void bindTexture (Shader* shader)
    {

        // Bind Textures using texture units
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture);

        glUniform1i(glGetUniformLocation(shader->Program, "material.diffuse"),  0);
        glUniform1i(glGetUniformLocation(shader->Program, "material.specular"), 1);
        glUniform1f(glGetUniformLocation(shader->Program, "material.shininess"), 32.0f);


    }

    void unbindTexture()
    {
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    GLuint& getTexture()
    {
        return texture;
    }


private:
    static unsigned int instanceID;
    GLuint texture;
    int wrapSParameter;
    int wrapTParameter;

    int minFilter;
    int magFilter;

    unsigned int id;

    int width, height;

    unsigned char* image;
    char* path;


    void loadImage()
    {
        glBindTexture(GL_TEXTURE_2D, texture);
        image = SOIL_load_image(path, &width, &height, 0, SOIL_LOAD_RGB);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
        glGenerateMipmap(GL_TEXTURE_2D);
        SOIL_free_image_data(image);



        glBindTexture(GL_TEXTURE_2D, 0);  // Unbind texture when done, so we won't accidentily mess up our texture.
    }


};


unsigned int MyTexture :: instanceID=0;
#endif
