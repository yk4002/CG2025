#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <glm/glm.hpp>
#include "Colour.h"
#include "ModelTriangle.h"

// Read MTL file and return a map from material name to Colour
std::unordered_map<std::string, Colour> readMtlFile(const std::string &filename);

// Read OBJ file and return a vector of ModelTriangles
std::vector<ModelTriangle> readObjFile(const std::string &objFilename, float scale);