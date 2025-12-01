#include "parser.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <Utils.h>

std::unordered_map<std::string, Colour> readMtlFile(const std::string &filename) {
    //use a hashmap or hashtable for more
    //efficient colour lookup (using the name of the colour as a key).
    //there are 3 numbers for each thing?
    std::unordered_map<std::string, Colour> hash;

    std::string textLine;
    std::ifstream MyReadFile(filename);
    while (getline(MyReadFile, textLine)) {

        //first check if newmtl is declared
        if (textLine.rfind("newmtl", 0) == 0)  {
            std::vector<std::string> array = split(textLine, ' ');
            std::string name = array[1]; //assign name

            //now somehow read the next line and convert each of those values into rgb
            if (getline(MyReadFile, textLine)) {
                //check for the KD that indicates colours
                if (textLine.rfind("Kd", 0) == 0) {
                    std::vector<std::string> c = split(textLine, ' ');
                    //multiply each thing by 255 to get its rgb value
                    int red = std::stof(c[1]) * 255;
                    int green = std::stof(c[2]) * 255;
                    int blue = std::stof(c[3]) * 255;
                    hash[name] = Colour(name, red, green, blue);
                }
            }
        }
    }

    return hash;
}




//use the data the mtl file,followed by the obj file contains to populate a vector of ModelTriangles
//maybe make it exclusively something that reads info
std::vector<ModelTriangle> readObjFile(const std::string &objFilename, float scale, const std::string &mtlFilename) {

    //read the mtl files and store in hashtable
    std::unordered_map<std::string, Colour> hash = readMtlFile(mtlFilename);

    //declare vertices and model triangle vectors
    std::vector<glm::vec3> vertices;
    std::vector<ModelTriangle> modTriVec;
    std::vector<TexturePoint> textPts;

    //declare colour variable
    Colour colour("default", 255, 255, 255);

    //readfile logic
    std::string textLine;
    std::ifstream MyReadFile(objFilename);
    while (getline(MyReadFile, textLine)) {

        // use material
        if (textLine.rfind("usemtl", 0) == 0) {
            auto splitVec = split(textLine, ' ');
            colour = hash[splitVec[1]];
        }

        // vertex position
        else if (textLine.rfind("v ", 0) == 0) {
            auto splitVec = split(textLine, ' ');
            float x = std::stof(splitVec[1]) * scale;
            float y = std::stof(splitVec[2]) * scale;
            float z = std::stof(splitVec[3]) * scale;
            vertices.emplace_back(x, y, z);
        }

        // texture coordinate
        else if (textLine.rfind("vt ", 0) == 0) {
            auto splitVec = split(textLine, ' ');
            float u = std::stof(splitVec[1]);
            float v = std::stof(splitVec[2]);
            textPts.emplace_back(TexturePoint(u, v));
        }


        // face loop
        else if (textLine.rfind("f ", 0) == 0) {
            auto splitVec = split(textLine, ' ');
            std::array<glm::vec3, 3> triVerts;
            std::vector<TexturePoint> triTexts;
            bool textureFace = false;


            //loop for each of the 3 trianglepoints
            for (int i = 0; i < 3; ++i) {
                //each v/vt is a token
                std::string vAndT = splitVec[i + 1];

                //split these based on slashes
                auto pair = split(vAndT, '/');

                //Read the vertex index and then access the value from vector
                int vIndex = std::stoi(pair[0]) - 1; //convert back to regular indexing
                triVerts[i] = vertices[vIndex];

                //Likewise for texture index if it exists
                if(pair.size() > 1 && !pair[1].empty()) {
                    int tIndex = std::stoi(pair[1]) - 1; //is this always viable
                    if (tIndex >= 0 && tIndex < textPts.size()) {
                        triTexts.push_back(textPts[tIndex]);
                    }

                }
                else {
                    //push negative texturepoint to indicate lack of texture
                    triTexts.push_back(TexturePoint{-1.0f, -1.0f});
                }
            }

            // build triangle after the loop
            ModelTriangle tri(triVerts[0], triVerts[1], triVerts[2], colour);
            tri.texturePoints = triTexts; //THIS ASSIGNMENT CAUSES SEGFAULT


            //calculate the triangle normal
            glm::vec3 v0 = tri.vertices[0];
            glm::vec3 v1 = tri.vertices[1];
            glm::vec3 v2 = tri.vertices[2];
            tri.normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));

            //add triangle
            modTriVec.push_back(tri);
        }
    }

    //close file and return modtriVec - CORRECT
    MyReadFile.close();
    return modTriVec;
}