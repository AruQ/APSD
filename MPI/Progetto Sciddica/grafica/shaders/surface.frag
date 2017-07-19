#version 330 core
struct Material {
    sampler2D diffuse;
    sampler2D specular;
    float     shininess;
};

struct Light {
    vec3 position;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;

    float constant;
    float linear;
    float quadratic;
};

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoords;

in float Red;

out vec4 color;

uniform vec3 viewPos;
uniform Material material;
uniform Light light;

void main()
{

    vec3 tex1/*=vec3(texture(material.diffuse, TexCoords))*/;
    vec3 tex2/* = vec3(texture(material.specular, TexCoords))*/;
    if (Red == 0.0f)
    {
        tex1 =  vec3(texture(material.diffuse, TexCoords));
        tex2 =  vec3(texture(material.specular, TexCoords));

    }
    else
    {
        tex1 = vec3(1.0f, Red-0.01f,0.0f);
        tex2 = vec3(1.0f, Red-0.01f,0.0f);
    }
    // Ambient
    vec3 ambient = light.ambient * tex1;

//    // Diffuse
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(light.position - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = light.diffuse * diff * tex1;

    // Specular
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), material.shininess);
    vec3 specular = light.specular * spec * tex2;

    // Attenuation
    float distance    = length(light.position - FragPos);
    float attenuation = 1.0f / (light.constant + light.linear * distance + light.quadratic * (distance * distance));

    ambient  *= attenuation;
    diffuse  *= attenuation;
    specular *= attenuation;

//    vec3 mah = vec3(texture(material.diffuse, TexCoords));
    color = vec4(ambient + diffuse + specular, 1.0f);
//    color = vec4(mah, 1.0f);
}


