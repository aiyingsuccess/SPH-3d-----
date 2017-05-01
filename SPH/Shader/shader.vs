void main()
{
	vec3 posEye = vec3(gl_ModelViewMatrix * vec4(gl_Vertex.xyz, 1.0));
    float dist = length(posEye);
	gl_PointSize = 250.0/dist;
   // gl_PointSize = 10;
	gl_Color=vec4(1.0f,0.0f,0.0f,1.0f);
	gl_TexCoord[0] = gl_MultiTexCoord0;
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
    gl_FrontColor = gl_Color;
}