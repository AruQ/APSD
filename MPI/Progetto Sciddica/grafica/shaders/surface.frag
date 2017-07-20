#version 330 core
struct Material {
    sampler2D diffuse;
    sampler2D specular;
    float     shininess;
};


in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoords;

in float Red;

out vec4 color;

uniform vec3 viewPos;
uniform Material material;


void main()
{

    vec3 tex1/*=vec3(texture(material.diffuse, TexCoords))*/;
   
    if (Red == 0.0f)
    {
        tex1 =  vec3(texture(material.diffuse, TexCoords));
       
    }
    else
    {
        tex1 = vec3(1.0f, Red-0.01f,0.0f);
       
    }
    
    
    float ambientStrength = 0.7f;
    vec3 ambient = ambientStrength * vec3(1.0f,1.0f,1.0f);
    vec3 result = ambient * tex1;
    color = vec4(result, 1.0f);

 
}


