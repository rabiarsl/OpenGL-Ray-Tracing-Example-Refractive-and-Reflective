// Homework4.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include "Angel.h"

using namespace std;

typedef Angel::vec3  color3;
typedef Angel::vec4  point4;
typedef Angel::vec3 point3;


struct Ray
{
    point3 center;
    vec3 direction;
};

struct Sphere {
    point3 center;
    float radius;
    int pigmentNum;
    int surfFinNum;
};

int width;
int height;
string output_file;
point3 camera;
point3 at;
point3 up;
float fovy;
float lightSources[20][9];
int numOfLights;

float pigments[25][3];
int numOfPigments;

float surfaceFinishes[25][7];
int numOfSurfFin;

Sphere objects[100];
int numOfObjects;

color3 background = color3(0.5, 0.5, 0.5);

Sphere currentSphere;
Sphere lastSphere;

bool isInside = false;

int max_depth = 4;

point3 intersect(Ray ray, bool *status) {
    point3 O = ray.center;
    vec3 d = ray.direction;
    float minT = INFINITY;
    for (int i = 0; i < numOfObjects; i++) {
        Sphere sphere = objects[i];
        point3 c = sphere.center;
        float r = sphere.radius;
        vec3 u = c - O;
        float A = dot(d, d);
        float B = -2 * dot(d, u);
        float C = dot(u, u) - r * r;
        float delta = B * B - 4 * A * C;
        float EPS = 0.001;
        float t = -1;
        int s = 0;
        
        if (-EPS < delta && delta < EPS) {
            if (-B / 2 * A > EPS) {
                t = -B / 2 * A;
                s = 1;
            }
        }
        else if (delta > EPS) {
            float t1 = (-B - sqrt(delta)) / (2 * A);
            float t2 = (-B + sqrt(delta)) / (2 * A);
            if (t1 <= EPS && t2 > EPS){
                t = t2;
                s = 2;
            }
            else if (t1 > EPS) {
                t = t1;
                s = 3;
            }
        }
        if (t > EPS) {
            *status = true;
            if (t < minT) {
                minT = t;
                currentSphere = sphere;
                if (s == 2)
                    isInside = true;
                else
                    isInside = false;
            }
        }
    }
    if (*status) {
        point3 p = O + minT * d;
        return p;
    }
    return NULL;
}

vec3 compute_normal(point3 p) {
    vec3 N;
    if (isInside)
        N = -normalize(p - currentSphere.center);
    else
        N = normalize(p - currentSphere.center);
    return N;
}

bool visible(point3 p, point3 light) {
    point3 O = light;
    vec3 d = normalize(p - light);
    float minT = INFINITY;
    float EPS = 0.01;
    for (int i = 0; i < numOfObjects; i++) {
        Sphere sphere = objects[i];
        point3 c = sphere.center;
        float r = sphere.radius;
        vec3 u = c - O;
        float A = dot(d, d);
        float B = -2 * dot(d, u);
        float C = dot(u, u) - r * r;
        float delta = B * B - 4 * A * C;
        float t = -1;
        if (-EPS < delta && delta < EPS) {
            if (-B / 2 * A > 0) {
                t = -B / 2 * A;
            }
        }
        else if (delta > EPS) {
            float t1 = (-B - sqrt(delta)) / (2 * A);
            float t2 = (-B + sqrt(delta)) / (2 * A);
            if (t1 <= EPS && t2 > EPS) {
                t = t2;
            }
            else if (t1 > EPS) {
                t = t1;
            }
        }
        if (t > EPS) {
            if (t < minT) {
                minT = t;
            }
        }
    }
    point3 p2 = O + minT * d;
    if (abs(p.x - p2.x) < EPS &&
        abs(p.y - p2.y) < EPS &&
        abs(p.z - p2.z) < EPS)
        return true;
    return false;

}

color3 phong(int i, point3 p, vec3 n) {
    point3 light = point3(lightSources[i][0], lightSources[i][1], lightSources[i][2]);
    vec3 l = normalize(light - p);
    vec3 v = normalize(camera - p);
    vec3 h = normalize(l + v);
    color3 Li = color3(lightSources[i][3], lightSources[i][4], lightSources[i][5]);
    float a = lightSources[i][6];
    float b = lightSources[i][7];
    float c = lightSources[i][8];
    int surf = currentSphere.surfFinNum;
    float ka = surfaceFinishes[surf][0];
    float kd = surfaceFinishes[surf][1];
    float ks = surfaceFinishes[surf][2];
    float alfa = surfaceFinishes[surf][3];
    float kr = surfaceFinishes[surf][4];
    int pigment = currentSphere.pigmentNum;
    color3 C = color3(pigments[pigment][0], pigments[pigment][1], pigments[pigment][2]);
    if (i == 0) {
        return ka * Li * C;
    }
    else {
        float d = sqrt((p.x - light.x) * (p.x - light.x) + (p.y - light.y) * (p.y - light.y) + (p.z - light.z) * (p.z - light.z));
        return (Li / (a + b * d + c * d * d)) * (kd * C * fmax(0.0, dot(n, l)) + ks * pow(fmax(0.0, dot(n, h)), alfa));
    }
}

Ray reflect(Ray r, point3 p, vec3 normal) {
    vec3 d = r.direction;
    vec3 v = -normalize(d);
    Ray ret;
    ret.center = p;
    ret.direction = normalize(2 * dot(normal, v) * normal - v);
    return ret;
}

Ray transmit(Ray r, point3 p, vec3 n, float n2) {
    Ray ret;
    ret.center = p;
    float n1 = 0.35;
    vec3 d = r.direction;
    vec3 v = -normalize(d);
    float cost1 = dot(v, n);
    vec3 u1 = dot(v, n) * n;
    float sint1 = sqrt(1 - cost1 * cost1);
    if (n1 / n2 * sint1 > 1)
        return reflect(r, p, n);
    vec3 w1 = u1 - v;
    vec3 w2 = (n1 / n2) * w1;
    float cost2 = sqrt(1 - (n1 / n2) * (n1 / n2) * (1 - dot(v, n) * dot(v, n)));
    vec3 t = ((n1 / n2) * cost1 - cost2) * n - (n1 / n2) * v;
    ret.direction = normalize(t);
    return ret;
}

color3 trace(Ray ray, int depth) {
    bool status = false;
    color3 localC, reflectedC, transmittedC;
    point3 p;
    vec3 normal;
    if (depth > max_depth)
        return background;
    p = intersect(ray, &status);
    if (!status)
        return background;
    normal = compute_normal(p);
    localC = color3(0, 0, 0);
    for (int i = 0; i < numOfLights; i++) {
        point3 light = point3(lightSources[i][0], lightSources[i][1], lightSources[i][2]);
        if (i==0 || visible(p, light)) {
            localC += phong(i, p, normal);
        }
    }
    int surf = currentSphere.surfFinNum;
    float kr = surfaceFinishes[surf][4];
    float kt = surfaceFinishes[surf][5];
    float n2 = surfaceFinishes[surf][6];
    if (kr > 0.001) {
        Ray Rr = reflect(ray, p, normal);
        reflectedC = trace(Rr, depth + 1);
    }
    if (kt > 0.001) {
        Ray Rt = transmit(ray, p, normal, n2);
        transmittedC = trace(Rt, depth + 1);
    }
    return localC + kr*reflectedC + kt*transmittedC;
}


int main(int argc, char** argv)
{
    string filename = "test2refractive.in";
    if (argc >= 2)
        filename = argv[1];
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) {
        cerr << "Unable to open file";
        exit(1);   // call system to stop
    }

    // READ OUTPUT NAME
    inFile >> output_file;;
    string temp = "";

    //WIDTH & HEIGHT
    inFile >> temp;
    width = atoi(temp.c_str());
    inFile >> temp;
    height = atoi(temp.c_str());

    // READ CAMERA
    inFile >> temp;
    camera.x = atof(temp.c_str());
    inFile >> temp;
    camera.y = atof(temp.c_str());
    inFile >> temp;
    camera.z = atof(temp.c_str());

    // READ AT
    inFile >> temp;
    at.x = atof(temp.c_str());
    inFile >> temp;
    at.y = atof(temp.c_str());
    inFile >> temp;
    at.z = atof(temp.c_str());
    inFile >> temp;

    // READ UP
    up.x = atof(temp.c_str());
    inFile >> temp;
    up.y = atof(temp.c_str());
    inFile >> temp;
    up.z = atof(temp.c_str());

    // READ FOVY
    inFile >> temp;
    fovy = atof(temp.c_str());

    // READ LIGHT SOURCES
    inFile >> temp;
    numOfLights = atoi(temp.c_str());
    for (int i = 0; i < numOfLights; i++) {
        for (int j = 0; j < 9; j++) {
            inFile >> temp;
            lightSources[i][j] = atof(temp.c_str());
        }
    }

    // READ PIGMENTS
    inFile >> temp;
    numOfPigments = atoi(temp.c_str());
    for (int i = 0; i < numOfPigments; i++) {
        inFile >> temp;
        if (strcmp(temp.c_str(), "solid") == 0) {
            for (int j = 0; j < 3; j++) {
                inFile >> temp;
                pigments[i][j] = atof(temp.c_str());
            }
        }
    }

    // READ SURFACE FINISHES
    inFile >> temp;
    numOfSurfFin = atoi(temp.c_str());
    for (int i = 0; i < numOfSurfFin; i++) {
        for (int j = 0; j < 7; j++) {
            inFile >> temp;
            surfaceFinishes[i][j] = atof(temp.c_str());
        }
    }

    // READ OBJECTS
    inFile >> temp;
    numOfObjects = atoi(temp.c_str());
    for (int i = 0; i < numOfObjects; i++) {
        Sphere sphere;
        for (int j = 0; j < 7; j++) {
            inFile >> temp;
            if (j == 0)
                sphere.pigmentNum = atoi(temp.c_str());
            else if (j == 1)
                sphere.surfFinNum = atoi(temp.c_str());
            if (j == 3)
                sphere.center.x = atoi(temp.c_str());
            if (j == 4)
                sphere.center.y = atoi(temp.c_str());
            if (j == 5)
                sphere.center.z = atoi(temp.c_str());
            if (j == 6)
                sphere.radius = atoi(temp.c_str());
        }
        objects[i] = sphere;
    }
    inFile.close();


    color3** colors = new color3 * [height];
    for (int i = 0; i < height; i++)
        colors[i] = new color3[width];
    point3 cO = camera;
    point3 cZ = -normalize(at - camera);
    point3 cX = normalize(cross(up, cZ));
    point3 cY = cross(cZ, cX);
    float h = 2 * tan(fovy / 2 * DegreesToRadians);
    float w = h * width / height;
    Ray ray;
    ray.center = camera;
    ofstream file;
    file.open(output_file, fstream::out | fstream::binary);
    file << "P6" << "\n";
    file << width << " ";
    file << height << "\n";
    file << "255" << "\n";

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            float px = w * j / width - w / 2;
            float py = -h * i / height + h / 2;
            point3 pix = cO + px * cX + py * cY - cZ;
            ray.direction = normalize(pix - ray.center);
            color3 color = trace(ray, 0);
            colors[i][j] = color;
            unsigned char c;
            float temp = (colors[i][j].x > 1) ? 1.0 : colors[i][j].x;
            temp = (temp < 0) ? 0.0 : temp;
            c = temp * 255.0;
            file << c;
            temp = (colors[i][j].y > 1) ? 1.0 : colors[i][j].y;
            temp = (temp < 0) ? 0.0 : temp;
            c = temp * 255.0;
            file << c;
            temp = (colors[i][j].z > 1) ? 1.0 : colors[i][j].z;
            temp = (temp < 0) ? 0.0 : temp;
            c = temp * 255.0;
            file << c;
        }
    }
    file.close();
}

