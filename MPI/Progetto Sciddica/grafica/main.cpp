// Std. Includes
#include <string>

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

#include<mpi.h>

// GLFW
#include <GLFW/glfw3.h>

// GL includes
#include "Shader.h"
#include "Camera.h"
#include "Surface.h"
#include "PointLight.h"

#include "Texture.h"
#include "Temperature.h"

// GLM Mathemtics
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// Other Libs
#include <SOIL.h>


// Properties
GLuint screenWidth = 1000, screenHeight = 1000;

// Function prototypes
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void Do_Movement();


//LightProperties
Intensity intensity = Intensity (glm::vec3 (0.7f, 0.7f, 0.7f),glm::vec3 (0.9f, 0.9f, 0.9f),glm::vec3(1.0f, 1.0f, 1.0f) );
PointLight pointLight = PointLight(1,glm::vec4(0.0f, 0.0f, -1.0f,1.0f),1.0f, 0.00014,0.0000007, intensity);

// Camera
Camera camera(glm::vec3(0.0f, 0.0f, -1.0f));

bool keys[1024];
GLfloat lastX = 400, lastY = 300;
bool firstMouse = true;

GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;

bool requestReceive = false;
int updated = 0;

int size;
double * temp;
MPI_Request req;


void receive (/*MPI_Comm & comm*/)
{
    MPI_Irecv(temp, size, MPI_DOUBLE, 0, 314, MPI_COMM_WORLD, &req);
    requestReceive = true;
    updated=0;
}


void runGraphics (/*MPI_Comm & comm*/)
{
        Surface altitude("../data/dem.asc");
//    Surface altitude("../data/prova.dat");

    altitude.initialize();
    size = altitude.getCoordinates().nCols * altitude.getCoordinates().nRows;
    temp= new double[size];

    camera.Front = glm::vec3(0.0f, -1.0f,-1.0f);
    // Init GLFW
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

    GLFWwindow* window = glfwCreateWindow(screenWidth, screenHeight, "ModelLoading - LearnOpenGL", NULL, NULL); // Windowed
    glfwMakeContextCurrent(window);

    // Set the required callback functions
    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // Options
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // Initialize GLEW to setup the OpenGL Function pointers
    glewExperimental = GL_TRUE;
    glewInit();

    // Define the viewport dimensions
    glViewport(0, 0, screenWidth, screenHeight);

    // Setup some OpenGL options
    glEnable(GL_DEPTH_TEST);

    // Setup and compile our shaders
    Shader shader("../grafica/shaders/surface.vs", "../grafica/shaders/surface.frag");


    //    Altitude altitude("../data/altitudes.dat");



    //    Altitude altitude("../data/DEM_test.dat");
    //        Altitude altitude("../data/DEM_Albano.asc");


    //    Matrix matrix("../data/lava.dat");

    //    altitude.addMatrix(matrix);
//    Temperature temperature ("../data/prova.dat");
        Temperature temperature ("../data/source.asc");

//    cout<<" dim della temperatura "<<temperature.nRows <<"   "<<temperature.nCols<<endl;


    altitude.setTemperature(&temperature);
    //        temperature.printColor();

    MyTexture texture("../grafica/images/texture.png");
    texture.setParameters(GL_REPEAT, GL_REPEAT,GL_NEAREST, GL_NEAREST);


    double* vertices = altitude.getVBOVertices();

    cout<<pointLight<<endl;


    unsigned int* indices = altitude.getEBO();

    GLuint VBO, VAO, EBO;

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    // Bind our Vertex Array Object first, then bind and set our buffers and pointers.
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, altitude.getSizeVBO()*sizeof(double), vertices, GL_STREAM_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, altitude.getSizeEBO()*sizeof(unsigned int), indices, GL_STATIC_DRAW);


    // Position attribute
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 9 * sizeof(double), (GLvoid*)0);
    glEnableVertexAttribArray(0);
    //normals
    glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, 9 * sizeof(double), (GLvoid*)(3 * sizeof(double)));
    glEnableVertexAttribArray(1);
    // TexCoord attribute
    glVertexAttribPointer(2, 2, GL_DOUBLE, GL_FALSE, 9 * sizeof(double), (GLvoid*)(6 * sizeof(double)));
    glEnableVertexAttribArray(2);
    // red value
    glVertexAttribPointer(3, 1, GL_DOUBLE, GL_FALSE, 9 * sizeof(double), (GLvoid*)(8 * sizeof(double)));
    glEnableVertexAttribArray(3);

    //                    if (temperature[globalIndex] != globalIndex)
    //                        cout<<temperature[globalIndex]<<" ";


    glBindVertexArray(0); // Unbind VAO





    int count = 0;

    // Game loop
    while(!glfwWindowShouldClose(window))
    {
        // Set frame time
        GLfloat currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // Check and call events
        glfwPollEvents();
        Do_Movement();

        // Clear the colorbuffer
        glClearColor(0.05f, 0.05f, 0.05f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        shader.Use();   // <-- Don't forget this one!

        //        pointLight.setPosition(camera.Position);
        pointLight.setPosition(camera.Position);

        pointLight.SetUniformData(&shader,"light");

        // Transformation matrices
        glm::mat4 projection = glm::perspective(camera.Zoom, (float)screenWidth/(float)screenHeight, 0.1f, 40000.0f);
        glm::mat4 view = camera.GetViewMatrix();
        glUniformMatrix4fv(glGetUniformLocation(shader.Program, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glUniformMatrix4fv(glGetUniformLocation(shader.Program, "view"), 1, GL_FALSE, glm::value_ptr(view));

        // Draw the loaded model
        glm::mat4 model ;

        //        model =  glm::scale(model,glm::vec3(0.5f,0.5f,0.5f));
        glUniformMatrix4fv(glGetUniformLocation(shader.Program, "model"), 1, GL_FALSE, glm::value_ptr(model));



        texture.bindTexture(&shader);
        glBindVertexArray(VAO);

        if (!requestReceive)
        {

            receive(/*comm*/);
            updated = 0;
            cout<<"aggiorno "<<endl;
        }

        if (updated == 0)
        {

            MPI_Status st;
            MPI_Test(&req,&updated,&st);
            if(updated == 1)
            {
                cout <<"test: "<< updated<<endl;
            }
        }

        if (updated==1)
        {
            cout<<"aggiorno il valori della temperatura e li passo al VBO"<<endl;
            altitude.updateTemperature(temp);


            vertices = altitude.getVBOVertices();

            glBufferData(GL_ARRAY_BUFFER, altitude.getSizeVBO()*sizeof(double), vertices, GL_STREAM_DRAW);

            updated = 0;
            requestReceive= false;

        }

        glDrawElements(GL_TRIANGLES, altitude.getSizeEBO(), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        texture.unbindTexture();

        // Swap the buffers
        glfwSwapBuffers(window);


    }


    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);

    glfwTerminate();

}


#pragma region "User input"

// Moves/alters the camera positions based on user input
void Do_Movement()
{
    // Camera controls
    if(keys[GLFW_KEY_W])
        camera.ProcessKeyboard(_FORWARD, deltaTime);
    if(keys[GLFW_KEY_S])
        camera.ProcessKeyboard(_BACKWARD, deltaTime);
    if(keys[GLFW_KEY_A])
        camera.ProcessKeyboard(_LEFT, deltaTime);
    if(keys[GLFW_KEY_D])
        camera.ProcessKeyboard(_RIGHT, deltaTime);
}

// Is called whenever a key is pressed/released via GLFW
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);

    if(action == GLFW_PRESS)
        keys[key] = true;
    else if(action == GLFW_RELEASE)
        keys[key] = false;
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if(firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    GLfloat xoffset = xpos - lastX;
    GLfloat yoffset = lastY - ypos;

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}

#pragma endregion
