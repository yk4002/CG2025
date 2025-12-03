//libraries that nead to be included
#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <CanvasPoint.h> //represents a point on the drawing canvas
#include <Colour.h>
#include <TextureMap.h>
#include <algorithm>
#include <ModelTriangle.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unistd.h>
#include <RayTriangleIntersection.h>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <csignal>
#include <chrono>
#include <thread>
//#include "parser/parser.h"


//------------------------------- -------------------------------------------------------------------------------------------------
//Parser


//read an MTL file and store colour name and rgb details (in colour constructor) in hashmap
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




//rewrite slightly
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
        bool sphereRead;

        if (textLine.rfind("o", 0) == 0) {
            auto splitVec = split(textLine, ' ');
            if (splitVec[1] == "sphere") 
            sphereRead = true;
            else sphereRead = false;
        }

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
                if (sphereRead){
                    triVerts[i].x +=0.6f;
                    triVerts[i].y +=0.0f;
                    triVerts[i].z +=-0.5f;                                     
                }

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
            tri.texturePoints = triTexts;


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


//global variables
glm::mat3 oriMat; //might put this into main loop eventually and adjust type signature of functions
glm::vec3 upSaved(0.0f, 1.0f, 0.0f); //used for rotation around x axis
float scale = 0.35;
//define the window
#define WIDTH 320
#define HEIGHT 240
//various names of the stuff
std::vector<ModelTriangle> triVec;
std::vector<ModelTriangle> oVec = readObjFile("cornell-box.obj", scale, "cornell-box.mtl");
std::vector<ModelTriangle> sVec = readObjFile("sphere.obj", scale, "sphere.mtl");
std::vector<ModelTriangle> tVec = readObjFile("textured-cornell-box.obj", scale, "textured-cornell-box.mtl");


float ambient = 0.1f;
glm::vec3 lightSource(0.0f, 0.8f, 0.0f); //position of the light source (or its center)



//all the variables for the conditional keypresses
//camera variables
static bool left = false;
static bool right = false;
static bool up = false;
static bool down = false;
static bool zoomIn = false; //q
static bool zoomOut = false; //e
static bool rightOrb = false; 
static bool downOrb = false; 
static bool leftOrb = false; 
static bool upOrb = false; 
static bool orbitVar = false; 
static bool complex = false; 
static bool upLight = false;
static bool downLight = false;
static bool leftLight = false;
static bool rightLight = false;
static bool outLight = false;
static bool inLight = false;

//switching between renders
static bool wireF = false; 
static bool rast = false;
static bool ray = false; 
static bool diffuse = true;
static bool hShad = false; 
static bool sShad = false; 
static bool gourad = false; 
static bool phong = false; 
static bool text = false; 
static bool sphere = false;
static bool mirror = false; 
//other
static bool running = false;
static bool recording = false;

//handling the SDL event depending on its type
void handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {

        //moving camera position
       if (event.key.keysym.sym == SDLK_a) {
       std::cout << "" << std::endl;
       left = true;
       }
       else if (event.key.keysym.sym == SDLK_d) {
       std::cout << "" << std::endl;
       right = true;
       }
       else if (event.key.keysym.sym == SDLK_w) {
       std::cout << "" << std::endl;
       up = true;
       }
       else if (event.key.keysym.sym == SDLK_s) {
       std::cout << "" << std::endl;
       down = true;
       }
       else if (event.key.keysym.sym == SDLK_q) {
       std::cout << "" << std::endl;
       zoomIn = true;
       }
       else if (event.key.keysym.sym == SDLK_e) {
       std::cout << "" << std::endl;
       zoomOut = true;
       }
       
        //moving light position
       if (event.key.keysym.sym == SDLK_j) {
       std::cout << "" << std::endl;
       leftLight = true;
       }
       else if (event.key.keysym.sym == SDLK_l) {
       std::cout << "" << std::endl;
       rightLight = true;
       }
       else if (event.key.keysym.sym == SDLK_i) {
       std::cout << "" << std::endl;
       upLight = true;
       }
       else if (event.key.keysym.sym == SDLK_k) {
       std::cout << "" << std::endl;
       downLight = true;
       }
       else if (event.key.keysym.sym == SDLK_8) {
       std::cout << "" << std::endl;
       inLight = true;
       }
       else if (event.key.keysym.sym == SDLK_9) {
       std::cout << "" << std::endl;
       outLight = true;
       }



        //orbiting keypresses
        else if (event.key.keysym.sym == SDLK_UP) {
            std::cout << "" << std::endl;
            upOrb = true;
        }
        else if (event.key.keysym.sym == SDLK_LEFT) {
            std::cout << "" << std::endl;
            leftOrb = true;
        }
        else if (event.key.keysym.sym == SDLK_RIGHT) {
            std::cout << "" << std::endl;
            rightOrb = true;
        }
       else if (event.key.keysym.sym == SDLK_DOWN) {
            std::cout << "" << std::endl;
            downOrb = true;
        }

        //now switching between different modes
       else if (event.key.keysym.sym == SDLK_b) {
            std::cout << "b" << std::endl;
            rast = false;
            ray = false;
            wireF = true;
        }
       else if (event.key.keysym.sym == SDLK_n) {
            std::cout << "s" << std::endl;
            ray = false;
            wireF = false;
            rast = true;

        }
       else if (event.key.keysym.sym == SDLK_m) {
            std::cout << "m" << std::endl;
            rast = false;
            ray = true;
            wireF = false;

        }

        //automatic orbit conditional on keypress to turn on and off
        else if (event.key.keysym.sym == SDLK_o) {
            if (orbitVar == false) {
                 std::cout << "orbit on" << std::endl;
                 complex = false;
                 orbitVar = true;
            }
            else {
                std::cout << "orbit off" << std::endl;
                orbitVar = false;
            }
        }

        //make sure this doesn't loop, and that in the main loop it sets itself to be false
        else if (event.key.keysym.sym == SDLK_p) {
            if (complex == false) {
                 std::cout << "animation looping" << std::endl;
                 complex = true;
            }
            else {
                std::cout << "Animation stopping" << std::endl;
                complex = false;
            }
        }


        //shadows, texture and reflection
        else if (event.key.keysym.sym == SDLK_z) {
            if (hShad == false) {
                 std::cout << "hard shadows on" << std::endl;
                 sShad = false;
                 hShad = true;
            }
            else {
                std::cout << "hard shadows off" << std::endl;
                hShad = false;
            }
        }
        
        
        else if (event.key.keysym.sym == SDLK_x) {
            if (sShad == false) {
                 std::cout << "soft shadows on" << std::endl;
                 hShad = false;
                 sShad = true;
            }
            else {
                std::cout << "soft shadows off" << std::endl;
                sShad = false;
        }    }    
        

        else if (event.key.keysym.sym == SDLK_r) {
            if (mirror == false) {
                 std::cout << "mirror on" << std::endl;
                 mirror = true;
            }
            else {
                std::cout << "mirror off" << std::endl;
                mirror = false;
            }
        }

        else if (event.key.keysym.sym == SDLK_1) ambient += 0.05;
        else if (event.key.keysym.sym == SDLK_2) ambient -= 0.05;

       else if (event.key.keysym.sym == SDLK_t) {
            if (text == false) {
            std::cout << "textured floor on" << std::endl;
            text = true;
            }
            else {
                std::cout << "textured floor off" << std::endl;
                text = false;
            }
        }

       else if (event.key.keysym.sym == SDLK_y) {
            if (sphere == false) {
            std::cout << "sphere added" << std::endl;
            sphere = true;
            }
            else {
                std::cout << "sphere removed" << std::endl;
                sphere = false;
            }
        }

        //lighting and shading
        else if (event.key.keysym.sym == SDLK_3) {
            if (diffuse == false) {
                 std::cout << "diffuse lighting on" << std::endl;
                 gourad = false;
                 phong = false;
                 diffuse = true;
            }
            else {
                std::cout << "diffuse lighting off" << std::endl;
                diffuse = false;
            }
        }

        else if (event.key.keysym.sym == SDLK_4) {
            if (gourad == false) {
                 std::cout << "gourad shading on" << std::endl;
                 diffuse = false;
                 phong = false;
                 gourad = true;
            }
            else {
                std::cout << "gourad shading off" << std::endl;
                gourad = false;
            }
        }

        else if (event.key.keysym.sym == SDLK_5) {
            if (phong == false) {
                 std::cout << "phong shading on" << std::endl;
                 diffuse = false;
                 gourad = false;
                 phong = true;
            }
            else {
                std::cout << "phong shading off" << std::endl;
                gourad = false;
            }
        }

        else if (event.key.keysym.sym == SDLK_ESCAPE) { // Press Esc to exit
            std::cout << "Exiting..." << std::endl;
            running = false;
        }

        //recording!
        else if (event.key.keysym.sym == SDLK_RETURN) {
            if (recording==false ) {
                recording = true;
                std::cout << "Recording started" << std::endl;
            }
            else recording = false;
            std::cout << "Recording ended" << std::endl;
        }




    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
       window.savePPM("output.ppm");
       window.saveBMP("output.bmp");
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------
//gives line pixels between two specified points
std::vector<CanvasPoint> giveLinePixels(const CanvasPoint& p1, const CanvasPoint& p2) {
    std::vector<CanvasPoint> linePixels;

    // distances between points
    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;
    float dz = p2.depth - p1.depth;

    //find the steps (maybe adjust slightly if needed)
    float steps = std::max(std::abs(dx), std::abs(dy));

    //interpolate using steps
    float xStep = dx / steps;
    float yStep = dy / steps;
    float zStep = dz / steps;

    //initialise values to start with
    float x = p1.x;
    float y = p1.y;
    float z = p1.depth;

    //for each step add a canvas point and then increment xyz by their step
    for (int i = 0; i <= std::round(steps); i++) {
        CanvasPoint c(std::round(x), std::round(y), z);
        linePixels.push_back(c);
        x += xStep;
        y += yStep;
        z += zStep;
    }
    return linePixels;
}




//given the start and end point of line on canvas, draw a line
void drawLine(DrawingWindow &window, CanvasPoint p1, CanvasPoint p2, Colour col, std::array<std::array<float, WIDTH>, HEIGHT> &depthBuffer) {
    std::vector<CanvasPoint> pixels = giveLinePixels(p1, p2); // give all the line pixels between p1 and p2
    uint32_t colour = (255 << 24) + (int(col.red) << 16) + (int(col.green) << 8) + int(col.blue); //pack the colour

    //for every pixel
    for (CanvasPoint& pixel : pixels) {

        //only plot the pixel if it is within bounds!
        if (pixel.x < 0 || pixel.x >= WIDTH || pixel.y < 0 || pixel.y >= HEIGHT) {continue;}
        //the initial fill when depth buffer is 0
        float depth = pixel.depth;
        if (depth > (depthBuffer[pixel.y][pixel.x])) {
            window.setPixelColour(std::round(pixel.x), std::round(pixel.y),colour);
            depthBuffer[pixel.y][pixel.x] = depth;
        }
    }
}



//Draw white outline triangle
//if we have to, pass model triangle into the thing as well to get the colour
void drawOutlineTriangle(DrawingWindow &window, CanvasTriangle tri, Colour col, std::array<std::array<float, WIDTH>, HEIGHT> &depthBuffer) {
    //extract vertices
    CanvasPoint v0 = tri.v0();
    CanvasPoint v1 = tri.v1();
    CanvasPoint v2 = tri.v2();
    //then draw white stroked triangle to outline the fill
    Colour white = Colour(255,255,255);
    drawLine(window, v0, v1, white, depthBuffer);
    drawLine(window, v1, v2, white, depthBuffer);
    drawLine(window, v2, v0, white, depthBuffer);
}




//rasterises triangle
void drawFillTriangle(DrawingWindow &window, CanvasTriangle tri, Colour col,
                      std::array<std::array<float, WIDTH>, HEIGHT> &depthBuffer) {

    // Sort vertices 
    CanvasPoint v0 = tri.v0();
    CanvasPoint v1 = tri.v1();
    CanvasPoint v2 = tri.v2();
    if (v1.y < v0.y) std::swap(v0, v1);
    if (v2.y < v0.y) std::swap(v0, v2);
    if (v2.y < v1.y) std::swap(v1, v2);

    int y0 = std::round(v0.y);
    int y1 = std::round(v1.y);
    int y2 = std::round(v2.y);


    //Compute intersection point for splitting flat-top/bottom
    float dy = v2.y - v0.y;
    float px = v0.x + ((v2.x - v0.x) * (v1.y - v0.y)) / dy;
    float pz = v0.depth + ((v2.depth - v0.depth) * (v1.y - v0.y)) / dy;

    //Top triangle
    float dy1 = v1.y - v0.y;
    if (dy1 <= 1) drawLine(window,CanvasPoint(v0.x, y0, v0.depth),CanvasPoint(px, y0, pz),col, depthBuffer);
    else{
        //define slopes and start value
        float invSlope1 = (v1.x - v0.x) / dy1;
        float invSlope2 = (px - v0.x) / dy1;
        float zSlope1 = (v1.depth - v0.depth) / dy1;
        float zSlope2 = (pz - v0.depth) / dy1;

        float xStart = v0.x;
        float xEnd   = v0.x;
        float zStart = v0.depth;
        float zEnd   = v0.depth;
        //loop over all y steps and draw a line between x start and xend, then update those values with the gradient increment             
        for (int y = y0; y < y1; y++) {
            drawLine(window,CanvasPoint(std::round(xStart), y, zStart),CanvasPoint(std::round(xEnd), y, zEnd),col, depthBuffer);
            xStart += invSlope1;
            xEnd   += invSlope2;
            zStart += zSlope1;
            zEnd   += zSlope2;
        }
    }
    //Bottom half triangle
    float dy2 = v2.y - v1.y;
    if (dy2 <= 1) drawLine(window,CanvasPoint(px, y1, pz),CanvasPoint(v2.x, y1, v2.depth),col, depthBuffer);
    else{
        float invSlope1 = (v2.x - v1.x) / dy2;
        float invSlope2 = (v2.x - px) / dy2;
        float zSlope1 = (v2.depth - v1.depth) / dy2;
        float zSlope2 = (v2.depth - pz) / dy2;

        float xStart = v2.x;
        float xEnd   = v2.x;
        float zStart = v2.depth;
        float zEnd   = v2.depth;

        for (int y = y2; y > y1; y--) {
            drawLine(window,
                     CanvasPoint(std::round(xStart), y, zStart),
                     CanvasPoint(std::round(xEnd), y, zEnd),
                     col, depthBuffer);
            xStart -= invSlope1;
            xEnd   -= invSlope2;
            zStart -= zSlope1;
            zEnd   -= zSlope2;
        }
    }

    //draw line at y=y1 along the middle from p to v1
    drawLine(window,CanvasPoint(px, y1, pz),CanvasPoint(v1.x, y1, v1.depth),col, depthBuffer);
}





//-------------------------------------------------------------------------------------------------------------------------
//WEEK 5

//translate camera pos given the camera pos and xyz translation
glm::vec3 translatePos(glm::vec3 camPos, float x, float y, float z) {
    glm::vec3 v (x, y, z); //create vector in world space
    //make it relative to current orientation
    // glm::vec3 translate = oriMat * v;
    return camPos + v;
}

//rotate a camera pos given the pos as well as the 3 columns of the rotation matrix
glm::vec3 rotateCamPos (glm::vec3 &camPos, glm::vec3 col1, glm::vec3 col2, glm::vec3 col3) {
    glm::mat3 rotate(col1, col2, col3);
    return rotate*camPos;
}


//finding orientation vectors
//do we need to normalise in every function?
glm::vec3 findFVec(glm::vec3 disVec) {
    return -glm::normalize(disVec); // normalise the vector from cam to obj and flip it for convention
}
glm::vec3 findRVec(glm::vec3 fVec){
    glm::vec3 u = glm::normalize(upSaved);
    return glm::normalize(glm::cross(u, fVec));
}
glm::vec3 findUVec(glm::vec3 fVec, glm::vec3 rVec){
    return glm::normalize(glm::cross(fVec, rVec));
}


//change orientation matrix
//we only do this once for each time orientation is changed!
void loadOriMat (glm::vec3 &camPos, glm::vec3 point) {
    glm::vec3 camToVertex = point - camPos;
    glm::vec3 fVec = findFVec(camToVertex);
    glm::vec3 rVec = findRVec(fVec);
    glm::vec3 uVec = findUVec(fVec, rVec);
    //load into orientation matrix
    oriMat = glm::mat3(rVec, uVec, fVec);
    upSaved = uVec; //update this
}



//Use a newly found orientation matrix to adjust a camera-to-vertex direction vectors
//when projecting vertices onto the image plane.
glm::vec3 adjustVertex(glm::vec3 &camPos, glm::vec3 &objPos) {
    //apply the new orientation matrix to all points to get all adjusted points
    glm::vec3 camToVertex = objPos - camPos;
    return camToVertex*oriMat;
}

//You should call your new lookAt function to keep camera focused on centre every time object is rotated
void lookAt(glm::vec3 &camPos) {
    //find and load the new orientation matrix now that the camPos has been changed.
    glm::vec3 centre(0.0f, 0.0f, 0.0f);
    loadOriMat (camPos, centre);
}

//---------------------------------------------------------------------------------------------------------------------------------




//---------------------------------------------------------------------------------------------------------------------------------
//projects 3d vertices onto image space, along with storing their depth
CanvasPoint projectVertexOntoCanvasPoint(glm::vec3 cameraPosition, float focalLength, glm::vec3 vertexPosition) {
    float x = vertexPosition.x;
    float y = vertexPosition.y;
    float depth = vertexPosition.z; //turn this into relative position?
    //convert to 2dpoints using the formula
    float u = std::round((focalLength * (x / -depth)) + (WIDTH / 2));
    float v = std::round(-(focalLength * (y / -depth)) + (HEIGHT / 2));
    //store depth as this
    return CanvasPoint(u, v, -1/depth);
}




//creates wireframe render of image
void wireFrameRender(const std::vector<ModelTriangle> &modTriVec, DrawingWindow &window, glm::vec3 camPos, float focalLength, std::array<std::array<float, WIDTH>, HEIGHT> &depthBuffer) {

    std::vector<ModelTriangle> copy = modTriVec;
    for (ModelTriangle &modTri : copy) {
        glm::vec3 v0 = adjustVertex(camPos, modTri.vertices[0]);
        glm::vec3 v1 = adjustVertex(camPos, modTri.vertices[1]);
        glm::vec3 v2 = adjustVertex(camPos, modTri.vertices[2]);

        //project each vertex onto the image plane
        CanvasPoint c0 = projectVertexOntoCanvasPoint(camPos, focalLength, v0);
        CanvasPoint c1 = projectVertexOntoCanvasPoint(camPos, focalLength, v1);
        CanvasPoint c2 = projectVertexOntoCanvasPoint(camPos, focalLength, v2);

        // Create a canvas triangle from the projected points
        CanvasTriangle tri{c0, c1, c2};

        // Draw filled triangle using the triangle's colour
        drawOutlineTriangle(window, tri, Colour(0,0,0), depthBuffer);
    }
}


//rasterises the image
void rasterRender(const std::vector<ModelTriangle> &modTriVec, DrawingWindow &window, glm::vec3 camPos, float focalLength, std::array<std::array<float, WIDTH>, HEIGHT> &depthBuffer) {

    std::vector<ModelTriangle> copy = modTriVec;
    for (ModelTriangle &modTri : copy) {
        glm::vec3 v0 = adjustVertex(camPos, modTri.vertices[0]);
        glm::vec3 v1 = adjustVertex(camPos, modTri.vertices[1]);
        glm::vec3 v2 = adjustVertex(camPos, modTri.vertices[2]);

        //project each vertex onto the image plane
        CanvasPoint c0 = projectVertexOntoCanvasPoint(camPos, focalLength, v0);
        CanvasPoint c1 = projectVertexOntoCanvasPoint(camPos, focalLength, v1);
        CanvasPoint c2 = projectVertexOntoCanvasPoint(camPos, focalLength, v2);

        // Create a canvas triangle from the projected points
        CanvasTriangle tri{c0, c1, c2};

        // Draw filled triangle using the triangle's colour
        drawFillTriangle(window, tri, modTri.colour, depthBuffer);
    }
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------
//RAY TRACING


RayTriangleIntersection getClosestValidIntersection(glm::vec3 &lightSourcePos, glm::vec3 &rayDirection, const std::vector<ModelTriangle> &triVec) {

    //these will be the params of the final RayTriangleIntersection object.
    glm::vec3 closestIntersection;
    float distanceFromCamera = std::numeric_limits<float>::infinity(); //initialise to infinity to start with
    ModelTriangle closestTriangle;
    size_t triIndex;


    //will search through the all of the triangles in the current scene and
    //return details of the closest intersected triangle
    size_t i = 0; //iterator for our loop
    bool intersected = false; //a variable to check whether we've actually had an intersection

    for (const ModelTriangle &triangle : triVec) {
        //matrix maths to get our possible solution
        //edges
        glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        glm::vec3 SPVector = lightSourcePos - triangle.vertices[0]; //s-p which is camera pos take away the point?
        glm::mat3 DEMatrix(-rayDirection, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        //possible solution is a vec3 that consists of:
    //  the absolute distance along the ray from the camera to the intersection point
        float t = possibleSolution[0];
    //    the proportional distance along the triangle's first edge that the intersection point occurs
        float u = possibleSolution[1];
    //   the proportional distance along the triangle's second edge that the intersection point occurs
        float v = possibleSolution[2];


        if ((t < distanceFromCamera) && (u >= 0.0f) && (u <= 1.0f) && (v >= 0.0f) && (v <= 1.0f) && (u + v <= 1.0f) && (t > 0.0f)) {
            closestIntersection = lightSourcePos + t * rayDirection;
            distanceFromCamera = t;
            closestTriangle = triangle;
            triIndex = i;
            intersected = true; //as soon as we hit a valid triangle make this variable true
        }
        i++;//increment the counter
    }

    //we only return the closest triangle (ie the smallest distance)
    if (intersected) {return RayTriangleIntersection(closestIntersection, distanceFromCamera, closestTriangle, triIndex);}

    //default constructor - the way we check it is to see if the distance from camera is equal to infinity
    else {
        glm::vec3 defaultVec(0.0f, 0.0f, 0.0f);
        return RayTriangleIntersection(defaultVec,std::numeric_limits<float>::infinity(),ModelTriangle(defaultVec, defaultVec, defaultVec, Colour(0,0,0)),std::numeric_limits<size_t>::max());
    }
}


//generates multi point light centred around this
//std::vector<glm::vec3> multiPointLight(glm::vec3 L, float r) {
//    int N = 10;
//
//    vol = N*N*N;
//    std::vector<glm::vec3> pts(vol);
//    for (int i = 0; i < N; i++)
//    for (int j = 0; j < N; j++)
//    for (int k = 0; k < N; k++) {
//        pts.push_back(L + p * r) //use indexing instead!
//    }
//}

//rewrite this in a way that makes sense
std::vector<glm::vec3> multiPointLight(glm::vec3 L, float r) {
    int N = 7;
    int noPoints = pow(N, 3);
    std::vector<glm::vec3> pts;
    pts.reserve(noPoints);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < 4; k++) {
                // Define each component of the point (X, Y, Z). Add small offset to make sure it is never zero
                float X = (i + 0.1f) / N * 2 - 1;
                float Y = (j + 0.1f) / N * 2 - 1;
                float Z = (k + 0.1f) / N * 2 - 1;

                // Create the point using X, Y, Z
                glm::vec3 p(X, Y, Z);

                // Check if the point is inside the unit sphere
                if (glm::length(p) <= 1.0f)
                    pts.push_back(L + p * r);  // Add the point to the list, scaled by radius r
            }
        }
    }

    return pts;
}






//for every vertex, look for the triangles that share that vertex
//average those face normals by doing vector sum then normalising
glm::vec3 findVertexNorm(const std::vector<ModelTriangle> &triVec, const glm::vec3 &vPos) {
    glm::vec3 vertexNormSum(0.0f);
    float diff = 1e-5f; //a tiny offset to make up for how you cant equate stuff?
    for (const ModelTriangle &triangle : triVec) {
        // check if any vertex of the triangle is basically equal to vpos
        if (glm::length(triangle.vertices[0] - vPos) < diff ||glm::length(triangle.vertices[1] - vPos) < diff ||glm::length(triangle.vertices[2] - vPos) < diff) {
            vertexNormSum += triangle.normal;
        }
    }
    return glm::normalize(vertexNormSum);
}


//checks to see if all texturepoints in the triangle are valid
bool checkTexture(const ModelTriangle &tri) {
    bool text = false;
    for (const TexturePoint &t : tri.texturePoints) {
        if (t.x < 0 || t.y < 0) text = false;
        else text= true; // all zero
    }
    return text;
}


//hard coded to blue rn - use this in main loop
void makeTriReflective(std::vector<ModelTriangle> &triVec, bool b) {
    Colour c(255,255,0);

        for (ModelTriangle &tri: triVec) {
            if (b == true) {
                if (tri.colour.red == c.red && tri.colour.green == c.green && tri.colour.blue == c.blue) tri.isMirror = true;
            }
            else tri.isMirror = false;
    }
}





//diffuse lighting
void rayTraceRenderr(const std::vector<ModelTriangle> &triVec,DrawingWindow &window,glm::vec3 &camPos,float &focalLength,glm::vec3 &lightSource)
 {
    float z = -focalLength; // z value of image plane
    uint32_t black = (255 << 24) + (int(0) << 16) + (int(0) << 8) + int(0);

    for (int u = 0; u < WIDTH; u++) {
        for (int v = 0; v < HEIGHT; v++) {
            // Convert pixel to 3D point on image plane
            float x = (u - WIDTH / 2) * (-z / focalLength);
            float y = -(v - HEIGHT / 2) * (-z / focalLength);
            glm::vec3 point3D(x, y, z);
            glm::vec3 rayDirection = glm::normalize(oriMat * point3D);

            // Find closest intersection with the scene
            RayTriangleIntersection closestIntersection = getClosestValidIntersection(camPos, rayDirection, triVec);

            if (closestIntersection.distanceFromCamera != std::numeric_limits<float>::infinity()) {
                glm::vec3 intersectPt = closestIntersection.intersectionPoint;
                ModelTriangle triangle = closestIntersection.intersectedTriangle;
                Colour colour = triangle.colour;
                glm::vec3 norm = glm::normalize(triangle.normal);

                float brightness = 1.0f;

                //TOGGLE diffuse lighting
                if (diffuse) {
                    glm::vec3 surfaceToLight = glm::normalize(lightSource - intersectPt);

                    // Proximity and angle of incidence lighting
                    float R = glm::length(lightSource - intersectPt);
                    float pi = 3.14159265f;
                    float prox = 12.0f / (4.0f * pi * R * R);
                    float aoi = glm::clamp(glm::dot(norm, surfaceToLight), 0.0f, 1.0f);
                    brightness = prox * aoi;

                    // Ambient lighting
                    if (brightness < ambient) brightness = ambient;
                }

                // Specular highlight (always calculated regardless of diffuse)
                glm::vec3 rayInc = glm::normalize(lightSource - intersectPt);
                glm::vec3 rayRefl = glm::normalize(rayInc - 2.0f * norm * glm::dot(rayInc, norm));
                glm::vec3 V = glm::normalize(camPos - intersectPt);
                float spec = glm::dot(rayRefl, V);
                float s = pow(glm::max(spec, 0.0f), 128.0f);

                // Final color
                int r = glm::clamp(int((colour.red  * brightness) + s * 250.0f), 0, 255);
                int g = glm::clamp(int((colour.green * brightness) + s * 250.0f), 0, 255);
                int b = glm::clamp(int((colour.blue  * brightness) + s * 250.0f), 0, 255);
                uint32_t c = (255 << 24) + (r << 16) + (g << 8) + b;
                window.setPixelColour(u, v, c);

                // Hard shadows
                if (hShad) {
                    glm::vec3 shadowRay = lightSource - intersectPt;
                    glm::vec3 shadowRayDir = glm::normalize(shadowRay);
                    glm::vec3 offsetIntersectPt = intersectPt + 0.001f * shadowRayDir;
                    float shadowRayLength = glm::length(shadowRay);
                    RayTriangleIntersection shadowIntersection = getClosestValidIntersection(offsetIntersectPt, shadowRayDir, triVec);
                    if (shadowIntersection.distanceFromCamera < shadowRayLength) window.setPixelColour(u, v, black);
                }
            }
        }
    }
}




//gouroud or phong aoi (comment out whichever one isn't being used)
void rayTraceRender(const std::vector<ModelTriangle> &triVec, DrawingWindow &window, glm::vec3 &camPos, float &focalLength, glm::vec3 &lightSource) {
    float z = -focalLength; // z value of image plane
    uint32_t black = (255 << 24) + (int(0) << 16) + (int(0) << 8) + int(0);

    // Load texture
    TextureMap textureImg = TextureMap("texture.ppm");

    //loop through image plane coords
    for (int u = 0; u < WIDTH; u++) {
        for (int v = 0; v < HEIGHT; v++) {
            // Convert pixel to 3D point on image plane
            float x = (u - WIDTH / 2) * (-z / focalLength);
            float y = -(v - HEIGHT / 2) * (-z / focalLength);
            glm::vec3 point3D(x, y, z);
            glm::vec3 rayDirection = glm::normalize(oriMat*point3D);
            //note that if we multiply other way round

            // Find closest intersection with the scene
            RayTriangleIntersection closestIntersection = getClosestValidIntersection(camPos, rayDirection, triVec);
            glm::vec3 intersectPt = closestIntersection.intersectionPoint;
            ModelTriangle triangle = closestIntersection.intersectedTriangle;
            Colour colour = triangle.colour;
            glm::vec3 surfaceToLight = glm::normalize(lightSource - intersectPt);
            glm::vec3 norm = glm::normalize(triangle.normal);
            

            // Soft shadow ray setup - turn this into a function conditional on keypress
            float shadWeight = 1.0f;
            std::vector<float> shadRayLengths;
            std::vector<RayTriangleIntersection> shadIntscts;
            int points;
            if (sShad) {
                float radius = 0.8f;
                std::vector<glm::vec3> multiSource = multiPointLight(lightSource, radius) ;
                points = multiSource.size();
                for (glm::vec3 &l : multiSource) {
                    glm::vec3 shadRay = l - intersectPt;
                    glm::vec3 shadRayDir = glm::normalize(shadRay);
                    glm::vec3 offsetIntersectPt = intersectPt + 0.001f * shadRayDir;
                    //these two variables are the most important
                    float shadRayLength = glm::length(shadRay);
                    RayTriangleIntersection shadIntsct = getClosestValidIntersection(offsetIntersectPt, shadRayDir, triVec);
                    shadRayLengths.push_back(shadRayLength);
                    shadIntscts.push_back(shadIntsct);
                }
            }


//           Only proceed if there was a valid intersection
            if (closestIntersection.distanceFromCamera != std::numeric_limits<float>::infinity()) {

                //loop through each light point to evaluate its shadow weight
                if (sShad) {
                    float shadHit = 0;
                    for(int i=0; i < points; i++) {
                        if (shadIntscts[i].distanceFromCamera < shadRayLengths[i]) shadHit ++;
                    }
                    float prop = shadHit/points;
                    shadWeight = 1.0 - prop;
                }

                // Compute barycentric coordinates
                glm::vec3 v0 = triangle.vertices[0];
                glm::vec3 v1 = triangle.vertices[1];
                glm::vec3 v2 = triangle.vertices[2];
                glm::vec3 bary = baryCoords(v0,v1,v2, intersectPt);
                float A = bary.x;
                float B = bary.y;
                float C = bary.z;
                
                // Set texture if enabled
                if (text && checkTexture(triangle)) {
                        glm::vec2 t0(triangle.texturePoints[0].x, triangle.texturePoints[0].y);
                        glm::vec2 t1(triangle.texturePoints[1].x, triangle.texturePoints[1].y);
                        glm::vec2 t2(triangle.texturePoints[2].x, triangle.texturePoints[2].y);
                        // Interpolate texture coordinates
                        glm::vec2 hitpoint = C*t0 + A*t1 + B*t2;
                        //calculate and clamp texture coordinates, and calculate index
                        int texX = hitpoint.x*((textureImg.width)-1);
                        int texY = hitpoint.y*((textureImg.height)-1);
                        texY = textureImg.height - 1 - texY;
                        texX = glm::clamp(texX, 0, int(textureImg.width) - 1);
                        texY = glm::clamp(texY, 0, int(textureImg.height) - 1);
                        int texIndex = texY * textureImg.width + texX;
                        // Extract color by shifting and adding 0s
                        uint32_t pixCol = textureImg.pixels[texIndex];
                        colour.red   = (pixCol >> 16) & 0xFF;
                        colour.green = (pixCol >> 8) & 0xFF;
                        colour.blue  = pixCol & 0xFF;
                }

//----------------------------------------------------------------------------------------------------------------
                //Mirror code, condense this again
                // Mirror code with soft shadows
if (mirror && triangle.isMirror) {
    glm::vec3 rayInc = rayDirection;
    glm::vec3 rayMirr = glm::normalize(rayInc - 2.0f * norm * glm::dot(rayInc, norm));
    glm::vec3 newOrigin = intersectPt + 0.01f * norm; // Slightly offset the point to avoid self-intersection

    // Trace the reflected ray fully
    RayTriangleIntersection i2 = getClosestValidIntersection(newOrigin, rayMirr, triVec);

    // Check if a valid intersection was found
    if (i2.distanceFromCamera != std::numeric_limits<float>::infinity()) {
        ModelTriangle tri2 = i2.intersectedTriangle;
        glm::vec3 hitPt2 = i2.intersectionPoint;

        // Compute barycentric coordinates for the reflected triangle
        glm::vec3 v0_2 = tri2.vertices[0];
        glm::vec3 v1_2 = tri2.vertices[1];
        glm::vec3 v2_2 = tri2.vertices[2];
        glm::vec3 bary2 = baryCoords(v0_2, v1_2, v2_2, hitPt2);
        float A2 = bary2.x;
        float B2 = bary2.y;
        float C2 = bary2.z;

        Colour col2 = tri2.colour;

        // Handle texture on reflected triangle (if present)
        if (checkTexture(tri2)) {
            glm::vec2 t0(tri2.texturePoints[0].x, tri2.texturePoints[0].y);
            glm::vec2 t1(tri2.texturePoints[1].x, tri2.texturePoints[1].y);
            glm::vec2 t2(tri2.texturePoints[2].x, tri2.texturePoints[2].y);
            glm::vec2 hitTex = C2 * t0 + A2 * t1 + B2 * t2;

            int texX = glm::clamp(int(hitTex.x * (textureImg.width - 1)), 0, int(textureImg.width) - 1);
            int texY = glm::clamp(int(hitTex.y * (textureImg.height - 1)), 0, int(textureImg.height) - 1);
            int texIndex = texY * textureImg.width + texX;

            uint32_t pixCol = textureImg.pixels[texIndex];
            col2.red   = (pixCol >> 16) & 0xFF;
            col2.green = (pixCol >> 8) & 0xFF;
            col2.blue  = pixCol & 0xFF;
        }

        // Compute vertex normals for the reflected triangle
        std::vector<glm::vec3> vertexNormals2(3);
        vertexNormals2[0] = findVertexNorm(triVec, v0_2);
        vertexNormals2[1] = findVertexNorm(triVec, v1_2);
        vertexNormals2[2] = findVertexNorm(triVec, v2_2);

        glm::vec3 interpNorm2 = glm::normalize(C2 * vertexNormals2[0] + A2 * vertexNormals2[1] + B2 * vertexNormals2[2]);
        glm::vec3 surfaceToLight2 = glm::normalize(lightSource - hitPt2);

        // Proximity factor (light falloff)
        float R2 = glm::length(lightSource - hitPt2);
        float pi = 3.14159265;
        float prox2 = 14.0f / (4.0f * pi * R2 * R2);

        // Diffuse reflection
        float aoi2 = glm::clamp(glm::dot(interpNorm2, surfaceToLight2), 0.0f, 1.0f);
        float brightness2 = aoi2 * prox2;

        // Ambient lighting
        if (brightness2 < ambient) brightness2 = ambient;

        // Specular reflection (Phong model)
        glm::vec3 rayRefl2 = glm::normalize(surfaceToLight2 - 2.0f * interpNorm2 * glm::dot(surfaceToLight2, interpNorm2));
        glm::vec3 V2 = glm::normalize(camPos - hitPt2);
        float spec2 = glm::dot(rayRefl2, V2);
        float s2 = glm::clamp(pow(spec2, 64.0f), 0.0f, 1.0f);

        // Final reflected color (combine diffuse and specular)
        int r2 = glm::clamp(int(col2.red * brightness2 + s2 * 255.0f), 0, 255);
        int g2 = glm::clamp(int(col2.green * brightness2 + s2 * 255.0f), 0, 255);
        int b2 = glm::clamp(int(col2.blue * brightness2 + s2 * 255.0f), 0, 255);
        colour = Colour(r2, g2, b2);

        // Handle soft shadows for the reflected ray
        if (sShad) {
            float shadWeight2 = 1.0f;
            std::vector<float> shadRayLengths2;
            std::vector<RayTriangleIntersection> shadIntscts2;
            int points2;

            float radius = 0.8f;
            std::vector<glm::vec3> multiSource2 = multiPointLight(lightSource, radius);
            points2 = multiSource2.size();

            // Loop through each light source for the reflected ray
            for (glm::vec3 &l : multiSource2) {
                glm::vec3 shadRay2 = l - hitPt2;
                glm::vec3 shadRayDir2 = glm::normalize(shadRay2);
                glm::vec3 offsetIntersectPt2 = hitPt2 + 0.001f * shadRayDir2;
                float shadRayLength2 = glm::length(shadRay2);
                RayTriangleIntersection shadIntsct2 = getClosestValidIntersection(offsetIntersectPt2, shadRayDir2, triVec);
                shadRayLengths2.push_back(shadRayLength2);
                shadIntscts2.push_back(shadIntsct2);
            }

            // Evaluate shadow hit proportion and adjust shadow weight
            float shadHit2 = 0;
            for (int i = 0; i < points2; i++) {
                if (shadIntscts2[i].distanceFromCamera < shadRayLengths2[i]) shadHit2++;
            }

            float prop2 = shadHit2 / points2;
            shadWeight2 = 1.0f - prop2;

            // Adjust final color based on shadow weight
            int rFinal = glm::clamp(int((colour.red * brightness2) * shadWeight2 + s2 * 255.0f), 0, 255);
            int gFinal = glm::clamp(int((colour.green * brightness2) * shadWeight2 + s2 * 255.0f), 0, 255);
            int bFinal = glm::clamp(int((colour.blue * brightness2) * shadWeight2 + s2 * 255.0f), 0, 255);

            // Set pixel color with soft shadow effect
            colour = Colour(rFinal, gFinal, bFinal);
        }
    }
    // Default to black if no valid intersection was found for reflection
    else colour = Colour(0, 0, 0);
}

//------------------------------------------------------------------------------------
                // Compute vertex normals
                std::vector<glm::vec3> vertexNormals(3);
                vertexNormals[0] = findVertexNorm(triVec, v0);
                vertexNormals[1] = findVertexNorm(triVec, v1);
                vertexNormals[2] = findVertexNorm(triVec, v2);
                glm::vec3 interpolatedNorm = glm::normalize(C * vertexNormals[0] + A * vertexNormals[1] + B * vertexNormals[2]);


                float aoi;
                //proximity
                float R = glm::length(lightSource - intersectPt);
                float pi = 3.14159265;
                float prox = 14.0f/(4.0f * pi * R * R);

                if (gourad) {
               //gourad shading
                   float br1 = glm::dot(vertexNormals[0], surfaceToLight);
                   float br2 = glm::dot(vertexNormals[1], surfaceToLight);
                   float br3 = glm::dot(vertexNormals[2], surfaceToLight);
                   float intBr = C*br1 + A*br2 + B*br3; //use vertex normals
                   aoi = glm::clamp(intBr, 0.0f, 1.0f);
                }
                else if (phong)  aoi = glm::clamp(glm::dot(interpolatedNorm, surfaceToLight), 0.0f, 1.0f);
                else if (diffuse) aoi = glm::clamp(glm::dot(norm, surfaceToLight), 0.0f, 1.0f);

                //create brightness variation as well as setting the ambient lighting
                float brightness = aoi*prox;
                //ambient lighting
                brightness = std::max(brightness, ambient);


                //specular lighting
                //if pink, then use interpolated normal from Phong
                glm::vec3 specN = norm;
                if (colour.red == 255.0f && colour.blue == 255.0f && colour.green == 0.0f && phong) specN = interpolatedNorm;
                glm::vec3 rayRefl = glm::normalize(surfaceToLight - 2.0f * specN * glm::dot(surfaceToLight, specN));
                glm::vec3 V = glm::normalize(camPos - intersectPt);
                float spec = glm::dot(rayRefl, V);
                float s = glm::clamp(pow(spec, 256.0f), 0.0f, 1.0f);


                //final colour
                int r = glm::clamp(int((colour.red  * brightness) * shadWeight + s*255.0f), 0, 255);
                int g = glm::clamp(int((colour.green * brightness) * shadWeight + s*255.0f), 0, 255);
                int b = glm::clamp(int((colour.blue  * brightness) * shadWeight + s*255.0f), 0, 255);
                uint32_t c = (255 << 24) + (r << 16) + (g << 8) + b;
                window.setPixelColour(u, v, c);

                // Hard shadows
                if (hShad) {
                    glm::vec3 shadowRay = lightSource - intersectPt;
                    glm::vec3 shadowRayDir = glm::normalize(shadowRay);
                    glm::vec3 offsetIntersectPt = intersectPt + 0.001f * shadowRayDir;
                    float shadowRayLength = glm::length(shadowRay);
                    RayTriangleIntersection shadowIntersection = getClosestValidIntersection(offsetIntersectPt, shadowRayDir, triVec);
                    if (shadowIntersection.distanceFromCamera < shadowRayLength) window.setPixelColour(u, v, black);
                }
            }
        }
    }
}











void record(DrawingWindow &window) {
    int frameCount = 0;

    // Create the output folder if needed (optional)
    system("mkdir -p frames");

    while (recording) {
        // Save each frame as a BMP image with a unique name
        std::ostringstream filename;
        filename << "frames/frame_" << std::setw(4) << std::setfill('0') << frameCount << ".bmp";
        window.saveBMP(filename.str().c_str());
        // Increment frame count
        frameCount++;
    }

    // After recording finishes, use FFmpeg to create the video from BMP files
    std::string command = "ffmpeg -framerate 30 -i frames/frame_%04d.bmp -c:v libx264 -r 30 -pix_fmt yuv420p output.mp4";
    system(command.c_str());

    std::cout << "Video saved as output.mp4" << std::endl;
}


//set triVec to be triggered by keypresses
void setTriVec() {
    if (sphere && text) {
        triVec = tVec;
        triVec.insert(triVec.end(), sVec.begin(), sVec.end());
    }
    else if (sphere) {
        triVec = oVec;
        triVec.insert(triVec.end(), sVec.begin(), sVec.end());  // Combine oVec and sVec
    }
    else if (text) triVec = tVec;
    else triVec = oVec;
}


//reset the depthBuffer to appropriate - check where this is called to try and solve the occlusion problem properly!
void depthBufReset(std::array<std::array<float, WIDTH>, HEIGHT> &depthBuffer) {
    for (int y = 0; y < HEIGHT; y++)
        for (int x = 0; x < WIDTH; x++)
            depthBuffer[y][x] = -std::numeric_limits<float>::infinity();
}




//general render function which changes based on the variables used
void render(const std::vector<ModelTriangle> &triVec, DrawingWindow &window, glm::vec3 &camPos, float focalLength, std::array<std::array<float, WIDTH>, HEIGHT> &depthBuffer, glm::vec3 lightSource) {
    if (wireF) {
        window.clearPixels();
        depthBufReset(depthBuffer);
        wireFrameRender(triVec, window, camPos, focalLength, depthBuffer);
    }
    else if (ray) {
        window.clearPixels();
        rayTraceRender(triVec, window, camPos, focalLength, lightSource);
    }
    else if (rast) {
        window.clearPixels();
        depthBufReset(depthBuffer);
        rasterRender(triVec, window, camPos, focalLength, depthBuffer);
}
}



//main function takes a pointer to argument array?
int main(int argc, char *argv[]) {
    srand(time(0));

    // Initialise important variables - adjustable
    glm::vec3 camPos(0.0f, 0.0f, 4.0f);  // Camera in front of the box
    float focalLength = 320.0f;          // Scale the image
    oriMat = glm::mat3(1.0f);  // Identity orientation matrix to begin with
    // depth buffer
    std::array<std::array<float, WIDTH>, HEIGHT> depthBuffer; 
    depthBufReset(depthBuffer);
    setTriVec(); //set the triangle vec to whatever is required and display an initial raster render


    //values related to translation and rotation
    float move = 0.1f;   // Translation step
    glm::vec3 centre(0.0f, 0.0f, 0.0f);
    float deg2rad = 3.14159265f / 180.0f;

    //get graphics working
    DrawingWindow window(WIDTH, HEIGHT, false);
    SDL_Event event;
    running = true;
    while (running) {
        if (window.pollForInputEvents(event)) {
            handleEvent(event, window);
        }
        rast = true; //initially set it to raster render
        render(triVec, window, camPos, focalLength, depthBuffer, lightSource);
        setTriVec();
        // record(window);

        //mirror toggling
         makeTriReflective(triVec, mirror);

         //phong sphere toggling

        // TRANSLATION
        if (left) {camPos = translatePos(camPos, -move, 0, 0);
            left = false;
        }
        if (right) {
            camPos = translatePos(camPos, move, 0, 0);
            right = false;
        }
        if (up) {
            camPos = translatePos(camPos, 0, move, 0);
            up = false;
        }
        if (down) {
            camPos = translatePos(camPos, 0, -move, 0);
            down = false;
        }
        if (zoomIn) {
            camPos = translatePos(camPos, 0, 0, -move);
            zoomIn = false;
        }
        if (zoomOut) {
            camPos = translatePos(camPos, 0, 0, move);
            zoomOut = false;
        }


        if (leftLight) {
            lightSource = translatePos(lightSource, -move, 0, 0);
            leftLight = false;
        }
        if (rightLight) {
            lightSource = translatePos(lightSource, move, 0, 0);
            rightLight = false;
        }
        if (upLight) {
            lightSource = translatePos(lightSource, 0, move, 0);
            upLight = false;
        }
        if (downLight) {
            lightSource = translatePos(lightSource, 0, -move, 0);
            downLight = false;
        }
        if (inLight) {
            lightSource = translatePos(lightSource, 0, 0, -move);
            inLight = false;
        }
        if (outLight) {
            lightSource = translatePos(lightSource, 0, 0, move);
            outLight = false;
        }

        // ROTATION
        if (leftOrb) {
            float theta = -5.0f * deg2rad;
            glm::vec3 c1(cos(theta), 0, -sin(theta));
            glm::vec3 c2(0, 1, 0);
            glm::vec3 c3(sin(theta), 0, cos(theta));
            camPos = rotateCamPos(camPos, c1, c2, c3);
            lookAt(camPos);
            leftOrb = false;
        }

        if (rightOrb) {
            float theta = 5.0f * deg2rad;
            glm::vec3 c1(cos(theta), 0, -sin(theta));
            glm::vec3 c2(0, 1, 0);
            glm::vec3 c3(sin(theta), 0, cos(theta));
            camPos = rotateCamPos(camPos, c1, c2, c3);
            lookAt(camPos);
            rightOrb = false;
        }

        if (upOrb) {
            float theta = 5.0f * deg2rad;
            glm::vec3 c1(1, 0, 0);
            glm::vec3 c2(0, cos(theta), -sin(theta));
            glm::vec3 c3(0, sin(theta), cos(theta));
            camPos = rotateCamPos(camPos, c1, c2, c3);
            lookAt(camPos);
            upOrb = false;
        }

        if (downOrb) {
            float theta = -5.0f * deg2rad;
            glm::vec3 c1(1, 0, 0);
            glm::vec3 c2(0, cos(theta), -sin(theta));
            glm::vec3 c3(0, sin(theta), cos(theta));
            camPos = rotateCamPos(camPos, c1, c2, c3);
            lookAt(camPos);
            downOrb = false;
        }

        // ORBIT FUNCTION (continuous Y rotation)
        if (orbitVar) {
            float theta = 0.5f * deg2rad;
            glm::vec3 c1(cos(theta), 0, -sin(theta));
            glm::vec3 c2(0, 1, 0);
            glm::vec3 c3(sin(theta), 0, cos(theta));
            camPos = rotateCamPos(camPos, c1, c2, c3);
            lookAt(camPos);
        }





        //Complex animation
    if (complex) {
    // Step 1: Move diagonally up (forward + upward)
    camPos = translatePos(camPos, 0.01f, 0.1f, 0.1f); // Move forward (1.0f in Z) and up (1.0f in Y)

    // // Step 2: Orbit around the scene
    // // Increment the rotation angle for a smooth rotation (let's rotate by 1 degree per frame)
    // float theta = 1.0f * deg2rad;

    // // Rotate the camera around the Y-axis (for horizontal orbit)
    // glm::vec3 c1(cos(theta), 0, -sin(theta)); // X-axis
    // glm::vec3 c2(0, 1, 0);  // Y-axis (center of rotation)
    // glm::vec3 c3(sin(theta), 0, cos(theta)); // Z-axis

    // // Apply the rotation
    // camPos = rotateCamPos(camPos, c1, c2, c3);

    // // Update the camera to look at the center
    // lookAt(camPos);

    // // End the animation after one full orbit (360 degrees)
    // static float orbitProgress = 0.0f;
    // orbitProgress += 5.0f; // Increment orbit progress (1 degree per frame)

    // if (orbitProgress >= 360.0f) {
    //     orbitProgress = 0.0f; // Reset orbit progress
    //     complex = false; // End the animation after a full 360-degree orbit
    // }
}



        // Render frame
        window.renderFrame();
    }
}
