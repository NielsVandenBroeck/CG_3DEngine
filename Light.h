//
// Created by Niels on 13/05/2021.
//

#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H

#include <vector>
#include "Libraries/vector3d.h"
using namespace std;
class Light {
public:
    Light(const vector<double> &ambientLight, const vector<double> &diffuseLight = {0.0,0.0,0.0}, const vector<double> &specularLight = {0.0,0.0,0.0})
            : ambientLight(ambientLight), diffuseLight(diffuseLight), specularLight(specularLight) {}
    vector<double> ambientLight;
    vector<double> diffuseLight;
    vector<double> specularLight;
    virtual Vector3D getlVector() {return Vector3D();}
    virtual Vector3D getLocation() {return Vector3D();}
    virtual double getAngle() {return 360;}
    virtual string getType() {return "ambient";}
};

class infLight: public  Light{
public:
    infLight(const vector<double> &ambientLight, const vector<double> &diffuseLight, const Vector3D &ldVector,
             const vector<double> &specularLight = {0.0,0.0,0.0}) : Light(ambientLight, diffuseLight,
                                                                                    specularLight),
                                                                              ldVector(ldVector) {}

    Vector3D getlVector() override {return ldVector;}
    string getType() override {return "inf";}

private:
    Vector3D ldVector;
};

class pointLight: public Light{
public:
    pointLight(const vector<double> &ambientLight, const vector<double> &diffuseLight,
               const Vector3D &location, double spotAngle,const vector<double> &specularLight  = {0.0,0.0,0.0}) : Light(ambientLight,
                                                                                                        diffuseLight,
                                                                                                        specularLight),
                                                                                                  location(location),
                                                                                                  spotAngle(
                                                                                                          spotAngle) {}
    Vector3D getLocation() override {return location;}
    double getAngle() override {return spotAngle;}
    string getType() override {return "point";}
private:
    Vector3D location;
    double spotAngle;
};

#endif //ENGINE_LIGHT_H
