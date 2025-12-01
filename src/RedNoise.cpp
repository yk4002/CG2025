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
//#include "parser/parser.h"

//define the window
#define WIDTH 320
#define HEIGHT 240


//global variables
glm::mat3 oriMat; //might put this into main loop eventually and adjust type signature of functions
glm::vec3 upSaved(0.0f, 1.0f, 0.0f); //used for rotation around x axis
std::vector<ModelTriangle> triVec;

//all the variables for the conditional keypresses
//camera variables
static bool left = false;
static bool right = false;
static bool up = false;
static bool down = false;
static bool zoomIn = false; //q
static bool zoomOut = false; //e
static bool rightOrb = false; //d
static bool downOrb = false; //s
static bool leftOrb = false; //a
static bool upOrb = false; //w
static bool orbitVar = false; //o
static bool complex = false; //p
//switching between renders
static bool wireF = false; //b
static bool rast = false; //n
static bool ray = false; //m
static bool hShad = false; //g
static bool sShad = false; //h
static bool gourad = false; //c
static bool phong = false; //v
static bool text = false; //j
static bool mirror = false; //k
static bool refr = false; //l
//other
static bool running = false;


//gives line pixels between two specified points
std::vector<CanvasPoint> giveLinePixels(const CanvasPoint& p1, const CanvasPoint& p2) {
    std::vector<CanvasPoint> line;

    // distances between points
    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;
    float dz = p2.depth - p1.depth;

    //find the steps (maybe adjust slightly if needed)
    float steps = std::max(std::abs(dx), std::abs(dy));

    //interpolate using steps
    float x_step = dx / steps;
    float y_step = dy / steps;
    float z_step = dz / steps;

    //initialise values to start with
    float x = p1.x;
    float y = p1.y;
    float z = p1.depth;

    //for each step add a canvas point and then increment xyz by their step
    for (int i = 0; i <= std::round(steps); i++) {
        CanvasPoint c((x), (y), z);
        line.push_back(c);
        x += x_step;
        y += y_step;
        z += z_step;
    }
    return line;
}




//given the start and end point of line on canvas, draw a line
void drawLine(DrawingWindow &window, CanvasPoint p1, CanvasPoint p2, Colour col, std::array<std::array<float, WIDTH>, HEIGHT> &depthBuffer) {
    std::vector<CanvasPoint> pixels = giveLinePixels(p1, p2); // give all the line pixels between p1 and p2
    uint32_t colour = (255 << 24) + (int(col.red) << 16) + (int(col.green) << 8) + int(col.blue); //pack the colour

    //for every pixel
    for (CanvasPoint& pixel : pixels) {

        //only plot the pixel if it is within bounds!
        if (pixel.x < 0 || pixel.x >= WIDTH || pixel.y < 0 || pixel.y >= HEIGHT) {continue;}
        float depth = pixel.depth;

        //the initial fill when depth buffer is 0
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
//fix the random horizontal lines
void drawFillTriangle(DrawingWindow &window, CanvasTriangle tri, Colour col,
                      std::array<std::array<float, WIDTH>, HEIGHT> &depthBuffer) {


    // Sort vertices by Y coordinate (v0 = top, v1 = middle, v2 = bottom)
    CanvasPoint v0 = tri.v0();
    CanvasPoint v1 = tri.v1();
    CanvasPoint v2 = tri.v2();
    if (v1.y < v0.y) std::swap(v0, v1);
    if (v2.y < v0.y) std::swap(v0, v2);
    if (v2.y < v1.y) std::swap(v1, v2);

    int y0 = std::round(v0.y);
    int y1 = std::round(v1.y);
    int y2 = std::round(v2.y);

    //height of 0? or 1?
    if (y0 == y2) {
        int x_min = std::round(std::min({v0.x, v1.x, v2.x}));
        int x_max = std::round(std::max({v0.x, v1.x, v2.x}));
        float z_min = std::min({v0.depth, v1.depth, v2.depth});
        float z_max = std::max({v0.depth, v1.depth, v2.depth});

        drawLine(window, CanvasPoint(x_min, y0, z_min),
                          CanvasPoint(x_max, y0, z_max),
                          col, depthBuffer);
        return;
    }

    //if height of triangle equals 1?
//    int height = y2-y1;
//    if (height = 1) {
//        window.setPixelColour();
//
//    }

    //or if width of triangle equals 1? Idk

    // Compute intersection point for splitting flat-top/bottom
    float dy_total = v2.y - v0.y;
    float px = v0.x + ((v2.x - v0.x) * (v1.y - v0.y)) / dy_total;
    float pz = v0.depth + ((v2.depth - v0.depth) * (v1.y - v0.y)) / dy_total;

    // Top flat-bottom triangle
    float dy1 = v1.y - v0.y;

    if (dy1 != 0) {
        float inv_slope_1 = (v1.x - v0.x) / dy1;
        float inv_slope_2 = (px - v0.x) / dy1;

        float z_slope_1 = (v1.depth - v0.depth) / dy1;
        float z_slope_2 = (pz - v0.depth) / dy1;

        //initialise these values
        float x_start = v0.x;
        float x_end   = v0.x;
        float z_start = v0.depth;
        float z_end   = v0.depth;

        for (int y = y0; y < y1; y++) {
            //why is this int suddenly?
            int xs = (x_start);
            int xe = (x_end);

            drawLine(window,
                     CanvasPoint(xs, y, z_start),
                     CanvasPoint(xe, y, z_end),
                     col, depthBuffer);

            x_start += inv_slope_1;
            x_end   += inv_slope_2;
            z_start += z_slope_1;
            z_end   += z_slope_2;
        }
    }



    //lower triangle
    float dy2 = v2.y - v1.y;

    if (dy2 != 0) {
        float inv_slope_1 = (v2.x - v1.x) / dy2;
        float inv_slope_2 = (v2.x - px) / dy2;

        float z_slope_1 = (v2.depth - v1.depth) / dy2;
        float z_slope_2 = (v2.depth - pz) / dy2;

        //initialise again, this time we are starting from the bottom
        float x_start = v2.x;
        float x_end   = v2.x;
        float z_start = v2.depth;
        float z_end   = v2.depth;

        for (int y = y2; y > y1; y--) {
            int xs = (x_start);
            int xe = (x_end);

            drawLine(window,
                     CanvasPoint(xs, y, z_start),
                     CanvasPoint(xe, y, z_end),
                     col, depthBuffer);

            x_start -= inv_slope_1;
            x_end   -= inv_slope_2;
            z_start -= z_slope_1;
            z_end   -= z_slope_2;
        }
    }


    //handle the case where y = y1 separately
    float x_left = px;
    float x_right = v1.x;
    float z_left = pz;
    float z_right = v1.depth;

    // Draw a single horizontal line at y1
    drawLine(window,
             CanvasPoint(x_left, y1, z_left),
             CanvasPoint(x_right, y1, z_right),
             col, depthBuffer);


}



//-------------------------------------------------------------------------------------------------------------------------
//WEEK 5

//translate camera pos given the camera pos and xyz translation
glm::vec3 translatePos(glm::vec3 camPos, float x, float y, float z) {
    glm::vec3 translate(x,y,z);
    return camPos + translate;
}

//rotate a camera pos given the pos as well as the 3 columns of the rotation matrix
glm::vec3 rotateCamPos (glm::vec3 &camPos, glm::vec3 col1, glm::vec3 col2, glm::vec3 col3) {
    glm::mat3 rotate(col1, col2, col3);
    return rotate*camPos;
}


//maybe add another function which spins the camera round in place
//does this need look at or nah? How is the orimat used?
//glm::vec3 spinCameraPos (glm::vec3 &camPos, glm::vec3 col1, glm::vec3 col2, glm::vec3 col3) {
//
//
//}

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
//put some texture code here

//------------------------------- -------------------------------------------------------------------------------------------------
//WEEK 4


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


//---------------------------------------------------------------------------------------------------------------------------------
//projects 3d vertices onto image space, along with storing their depth
CanvasPoint projectVertexOntoCanvasPoint(glm::vec3 cameraPosition, float focalLength, glm::vec3 vertexPosition) {
    float x = vertexPosition.x;
    float y = vertexPosition.y;
    float depth = vertexPosition.z; //turn this into relative position?
    //convert to 2dpoints using the formula
    float u = (focalLength * (x / -depth)) + (WIDTH / 2);
    float v = -(focalLength * (y / -depth)) + (HEIGHT / 2);
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
    int N = 10;
    std::vector<glm::vec3> pts;
    pts.reserve(N*N*N);

    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++) {
        glm::vec3 p(
            (i + 0.5f) / N * 2 - 1,
            (j + 0.5f) / N * 2 - 1,
            (k + 0.5f) / N * 2 - 1
        );
        if (glm::length(p) <= 1.0f)
            pts.push_back(L + p * r);
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
    Colour c(0,0,255);

        for (ModelTriangle &tri: triVec) {
            if (b == true) {
                if (tri.colour.red == c.red && tri.colour.green == c.green && tri.colour.blue == c.blue) tri.isMirror = true;
            }
            else tri.isMirror = false;
    }
}

//do something similar for refraction?




//diffuse lighting
void rayTraceRenderr(const std::vector<ModelTriangle> &triVec, DrawingWindow &window, glm::vec3 &camPos, float &focalLength, glm::vec3 &lightSource) {
    float z = -focalLength; // z value of image plane
    uint32_t black = (255 << 24) + (int(0) << 16) + (int(0) << 8) + int(0);

    for (int u = 0; u < WIDTH; u++) {
        for (int v = 0; v < HEIGHT; v++) {
            // Convert pixel to 3D point on image plane
            float x = (u - WIDTH / 2) * (-z / focalLength);
            float y = -(v - HEIGHT / 2) * (-z / focalLength);
            glm::vec3 point3D(x, y, z);
            glm::vec3 rayDirection = glm::normalize(oriMat*point3D);

            // Find closestintersection with the sscene
            RayTriangleIntersection closestIntersection = getClosestValidIntersection(camPos, rayDirection, triVec);
            glm::vec3 intersectPt = closestIntersection.intersectionPoint;
            ModelTriangle triangle = closestIntersection.intersectedTriangle;
            Colour colour = triangle.colour;
            glm::vec3 surfaceToLight = glm::normalize(lightSource - intersectPt);
            glm::vec3 norm = glm::normalize(triangle.normal);


            // Hard shadow ray setup
            glm::vec3 shadowRay = lightSource - intersectPt;
            glm::vec3 shadowRayDir = glm::normalize(shadowRay);
            glm::vec3 offsetIntersectPt = intersectPt + 0.001f * shadowRayDir;
            float shadowRayLength = glm::length(shadowRay);
            RayTriangleIntersection shadowIntersection = getClosestValidIntersection(offsetIntersectPt, shadowRayDir, triVec);

            //only if ray intersection exists
            if (closestIntersection.distanceFromCamera != std::numeric_limits<float>::infinity()) {
                float shadWeight = 1.0f;

                //proximity and angle of incidence lighting
                float R = glm::length(lightSource - intersectPt);
                float pi = 3.14159265;
                float prox = 12/(4.0f * pi * R * R);
                float aoi = glm::clamp(glm::dot(norm, surfaceToLight), 0.0f, 1.0f);
                float brightness = prox*aoi;
                //ambient lighting
                float ambient = 0.1f;
                if (brightness < ambient) brightness = ambient;


                //adding specular highlight
                glm::vec3 rayInc = surfaceToLight;
                glm::vec3 rayRefl = glm::normalize(rayInc - 2.0f * norm * glm::dot(rayInc, norm));
                glm::vec3 V = glm::normalize(camPos - intersectPt);
                float spec = glm::dot(rayRefl, V);
                float s = pow(spec, 128.0f); //raise it to the power of something to make area smaller

                //final colour
                int r = glm::clamp(int((colour.red  * brightness) + s*250.0f), 0, 255);
                int g = glm::clamp(int((colour.green * brightness) + s*250.0f), 0, 255);
                int b = glm::clamp(int((colour.blue  * brightness) + s*250.0f), 0, 255);
                uint32_t c = (255 << 24) + (r << 16) + (g << 8) + b;
                window.setPixelColour(u, v, c);

                //add hard shadows
                if (shadowIntersection.distanceFromCamera < shadowRayLength) {
                    window.setPixelColour(u,v,black);
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
//            float radius = 0.8f;
//            std::vector<glm::vec3> multiSource = multiPointLight(lightSource, radius) ;
//            int points = multiSource.size();
//            std::vector<float> shadRayLengths;
//            std::vector<RayTriangleIntersection> shadIntscts;
//            for (glm::vec3 &l : multiSource) {
//                glm::vec3 shadRay = l - intersectPt;
//                glm::vec3 shadRayDir = glm::normalize(shadRay);
//                glm::vec3 offsetIntersectPt = intersectPt + 0.001f * shadRayDir;
//                //these two variables are the most important
//                float shadRayLength = glm::length(shadRay);
//                RayTriangleIntersection shadIntsct = getClosestValidIntersection(offsetIntersectPt, shadRayDir, triVec);
//                shadRayLengths.push_back(shadRayLength);
//                shadIntscts.push_back(shadIntsct);
//            }


            // Only proceed if there was a valid intersection
//            if (closestIntersection.distanceFromCamera != std::numeric_limits<float>::infinity()) {
//
//                //loop through each light point to evaluate its shadow weight
//                float shadHit = 0;
//                for(int i=0; i < points; i++) {
//                    if (shadIntscts[i].distanceFromCamera < shadRayLengths[i]) shadHit ++;
//                }
//                float prop = shadHit/points;
//                float shadWeight = 1.0 - prop;


                // Compute barycentric coordinates
                glm::vec3 v0 = triangle.vertices[0];
                glm::vec3 v1 = triangle.vertices[1];
                glm::vec3 v2 = triangle.vertices[2];
                glm::vec3 bary = baryCoords(v0,v1,v2, intersectPt);
                float A = bary.x;
                float B = bary.y;
                float C = bary.z;

                 // Set color depending on texture
                if (checkTexture(triangle)) {
                    glm::vec2 t0(triangle.texturePoints[0].x, triangle.texturePoints[0].y);
                    glm::vec2 t1(triangle.texturePoints[1].x, triangle.texturePoints[1].y);
                    glm::vec2 t2(triangle.texturePoints[2].x, triangle.texturePoints[2].y);
                    // Interpolate texture coordinates
                    glm::vec2 hitpoint = C*t0 + A*t1 + B*t2;
                    //calculate and clamp texture coordinates, and calculate index
                    int texX = hitpoint.x*((textureImg.width)-1);
                    int texY = hitpoint.y*((textureImg.height)-1);
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
                if (triangle.isMirror) {
                    glm::vec3 rayInc = rayDirection;
                    glm::vec3 rayMirr = glm::normalize(rayInc - 2.0f * norm * glm::dot(rayInc, norm));
                    glm::vec3 newOrigin = intersectPt + 0.01f * norm;
                    // Trace the reflected ray fully
                    RayTriangleIntersection i2 = getClosestValidIntersection(newOrigin, rayMirr, triVec);

                    //loop 2 of ray trace
                    if (i2.distanceFromCamera != std::numeric_limits<float>::infinity()) {
                        ModelTriangle tri2 = i2.intersectedTriangle;
                        glm::vec3 hitPt2 = i2.intersectionPoint;

                        // compute barycentric coordinates
                        glm::vec3 v0_2 = tri2.vertices[0];
                        glm::vec3 v1_2 = tri2.vertices[1];
                        glm::vec3 v2_2 = tri2.vertices[2];
                        glm::vec3 bary2 = baryCoords(v0_2, v1_2, v2_2, hitPt2);
                        float A2 = bary2.x;
                        float B2 = bary2.y;
                        float C2 = bary2.z;

                        Colour col2 = tri2.colour;

                        // handle texture on reflected triangle
                        if (checkTexture(tri2)) {
                            glm::vec2 t0(tri2.texturePoints[0].x, tri2.texturePoints[0].y);
                            glm::vec2 t1(tri2.texturePoints[1].x, tri2.texturePoints[1].y);
                            glm::vec2 t2(tri2.texturePoints[2].x, tri2.texturePoints[2].y);
                            glm::vec2 hitTex = C2*t0 + A2*t1 + B2*t2;

                            int texX = glm::clamp(int(hitTex.x * (textureImg.width - 1)), 0, int(textureImg.width) - 1);
                            int texY = glm::clamp(int(hitTex.y * (textureImg.height - 1)), 0, int(textureImg.height) - 1);
                            int texIndex = texY * textureImg.width + texX;

                            uint32_t pixCol = textureImg.pixels[texIndex];
                            col2.red   = (pixCol >> 16) & 0xFF;
                            col2.green = (pixCol >> 8) & 0xFF;
                            col2.blue  = pixCol & 0xFF;
                        }

                        // compute vertex normals
                        std::vector<glm::vec3> vertexNormals2(3);
                        vertexNormals2[0] = findVertexNorm(triVec, v0_2);
                        vertexNormals2[1] = findVertexNorm(triVec, v1_2);
                        vertexNormals2[2] = findVertexNorm(triVec, v2_2);

                        glm::vec3 interpNorm2 = glm::normalize(C2 * vertexNormals2[0] + A2 * vertexNormals2[1] + B2 * vertexNormals2[2]);
                        glm::vec3 surfaceToLight2 = glm::normalize(lightSource - hitPt2);

                        //prox
                        float R2 = glm::length(lightSource - hitPt2);
                        float pi = 3.14159265;
                        float prox2 = 14.0f/(4.0f * pi * R2 * R2);
                        // diffuse
                        float aoi2 = glm::clamp(glm::dot(interpNorm2, surfaceToLight2), 0.0f, 1.0f);
                        float brightness2 = aoi2 * prox2;
                        //ambient
                        float ambient2 = 0.1f;
                        if (brightness2 < ambient2) brightness2 = ambient2;
                        // specular
                        glm::vec3 rayRefl2 = glm::normalize(surfaceToLight2 - 2.0f * interpNorm2 * glm::dot(surfaceToLight2, interpNorm2));
                        glm::vec3 V2 = glm::normalize(camPos - hitPt2);
                        float spec2 = glm::dot(rayRefl2, V2);
                        float s2 = glm::clamp(pow(spec2, 64.0f), 0.0f, 1.0f);

                        // final reflected color is calculated
                        int r2 = glm::clamp(int(col2.red   * brightness2 + s2 * 255.0f), 0, 255);
                        int g2 = glm::clamp(int(col2.green * brightness2 + s2 * 255.0f), 0, 255);
                        int b2 = glm::clamp(int(col2.blue  * brightness2 + s2 * 255.0f), 0, 255);
                        colour = Colour(r2, g2, b2);
                    }
                    //default to black
                    else colour = Colour(0, 0, 0);
                }
//------------------------------------------------------------------------------------
                // Compute vertex normals
                std::vector<glm::vec3> vertexNormals(3);
                vertexNormals[0] = findVertexNorm(triVec, v0);
                vertexNormals[1] = findVertexNorm(triVec, v1);
                vertexNormals[2] = findVertexNorm(triVec, v2);

                float aoi;
                //proximity
                float R = glm::length(lightSource - intersectPt);
                float pi = 3.14159265;
                float prox = 14.0f/(4.0f * pi * R * R);

//                //gourad shading
//                    float br1 = glm::dot(vertexNormals[1], surfaceToLight);
//                    float br2 = glm::dot(vertexNormals[2], surfaceToLight);
//                    float br3 = glm::dot(vertexNormals[3], surfaceToLight);
//                    float intBr = C*br1 + A*br2 + B*br3; //use vertex normals
//                    aoi = glm::clamp(intBr, 0.0f, 1.0f);

                // Phong shading
                glm::vec3 interpolatedNorm = glm::normalize(C * vertexNormals[0] + A * vertexNormals[1] + B * vertexNormals[2]);
                aoi = glm::clamp(glm::dot(interpolatedNorm, surfaceToLight), 0.0f, 1.0f);
                float brightness = aoi*prox;

                //ambient lighting
                float ambient = 0.2f;
                if (brightness < ambient) brightness = ambient;

                //specular lighting - not too sure if this is working please check
                //consider on how you add this term to the rest of the stuff
                glm::vec3 rayRefl = glm::normalize(surfaceToLight - 2.0f * norm * glm::dot(surfaceToLight, norm));
                glm::vec3 V = glm::normalize(camPos - intersectPt);
                float spec = glm::dot(rayRefl, V);
                float s = glm::clamp(pow(spec, 128.0f), 0.0f, 1.0f);
                s=0.0f;

                //this will eventually be used
                Colour finalColour;
                //final colour
                int r = glm::clamp(int(colour.red  * brightness + s*255.0f), 0, 255);
                int g = glm::clamp(int(colour.green * brightness + s*255.0f), 0, 255);
                int b = glm::clamp(int(colour.blue  * brightness + s*255.0f), 0, 255);
                uint32_t c = (255 << 24) + (r << 16) + (g << 8) + b;
                window.setPixelColour(u, v, c);


        }
    }
}





////multi point soft shadows
//void rayTraceRendesr(const std::vector<ModelTriangle> &triVec, DrawingWindow &window, glm::vec3 &camPos, float &focalLength, glm::vec3 &lightSource) {
//    float z = -focalLength; // z value of image plane
//    uint32_t black = (255 << 24) + (int(0) << 16) + (int(0) << 8) + int(0);
//
//    for (int u = 0; u < WIDTH; u++) {
//        for (int v = 0; v < HEIGHT; v++) {
//            // Convert pixel to 3D point on image plane
//            float x = (u - WIDTH / 2) * (-z / focalLength);
//            float y = -(v - HEIGHT / 2) * (-z / focalLength);
//            glm::vec3 point3D(x, y, z);
//            glm::vec3 rayDirection = glm::normalize(point3D);
//
//            // Find closest intersection with the scene
//            RayTriangleIntersection closestIntersection = getClosestValidIntersection(camPos, rayDirection, triVec);
//            glm::vec3 intersectPt = closestIntersection.intersectionPoint;
//            ModelTriangle triangle = closestIntersection.intersectedTriangle;
//            Colour colour = triangle.colour;
//            glm::vec3 surfaceToLight = glm::normalize(lightSource - intersectPt);
//            glm::vec3 norm = glm::normalize(triangle.normal);
//
//            // Soft shadow ray setup
//            float radius = 0.8f;
//            std::vector<glm::vec3> multiSource = multiPointLight(lightSource, radius) ;
//            int points = multiSource.size();
//            std::vector<float> shadRayLengths;
//            std::vector<RayTriangleIntersection> shadIntscts;
//            for (glm::vec3 &l : multiSource) {
//                glm::vec3 shadRay = l - intersectPt;
//                glm::vec3 shadRayDir = glm::normalize(shadRay);
//                glm::vec3 offsetIntersectPt = intersectPt + 0.001f * shadRayDir;
//                //these two variables are the most important
//                float shadRayLength = glm::length(shadRay);
//                RayTriangleIntersection shadIntsct = getClosestValidIntersection(offsetIntersectPt, shadRayDir, triVec);
//                shadRayLengths.push_back(shadRayLength);
//                shadIntscts.push_back(shadIntsct);
//            }
//
//
//            // Only proceed if there was a valid intersection
//            if (closestIntersection.distanceFromCamera != std::numeric_limits<float>::infinity()) {
//
//                //loop through each light point to evaluate its shadow weight
//                float shadHit = 0;
//                for(int i=0; i < points; i++) {
//                    if (shadIntscts[i].distanceFromCamera < shadRayLengths[i]) shadHit ++;
//                }
//                float prop = shadHit/points;
//                shadWeight = 1.0 - prop;
//
//
//                // Compute barycentric coordinates
//                glm::vec3 v0 = triangle.vertices[0];
//                glm::vec3 v1 = triangle.vertices[1];
//                glm::vec3 v2 = triangle.vertices[2];
//                glm::vec3 bary = baryCoords(v0,v1,v2, intersectPt);
//                float A = bary.x;
//                float B = bary.y;
//                float C = bary.z;
//
//                // Compute vertex normals
//                std::vector<glm::vec3> vertexNormals(3);
//                vertexNormals[0] = findVertexNorm(triVec, v0);
//                vertexNormals[1] = findVertexNorm(triVec, v1);
//                vertexNormals[2] = findVertexNorm(triVec, v2);
//
//                float brightness;
//
//                //gourad shading
//
////                    float br1aoi = glm::dot(vertexNormals[0], surfaceToLight);
////                    float br2aoi = glm::dot(vertexNormals[1], surfaceToLight);
////                    float br3aoi = glm::dot(vertexNormals[2], surfaceToLight);
////                    float aoi = C*br1 + A*br2 + B*br3; //use vertex normals
////                    //think it is abc but not sure
////                    brightness = glm::clamp(aoi, 0.0f, 1.0f);
//
//                // Phong shading with aoi
////                glm::vec3 interpolatedNorm = glm::normalize(C * vertexNormals[0] + A * vertexNormals[1] + B * vertexNormals[2]);
////                brightness = glm::clamp(glm::dot(interpolatedNorm, surfaceToLight), 0.0f, 1.0f);
//
//                //ambient lighting
//                float ambient = 0.1f;
//                if (brightness < ambient) brightness = ambient;
//
//                //specular lighting
//                glm::vec3 rayRefl = glm::normalize(surfaceToLight - 2.0f * norm * glm::dot(surfaceToLight, norm));
//                glm::vec3 V = glm::normalize(camPos - intersectPt);
//                float s = glm::clamp(pow(glm::dot(rayRefl, V), 64.0f), 0.0f, 1.0f);
//
//                //final colour (with soft shadow influences
//                int r = glm::clamp(int((colour.red  * brightness + s*255.0f) * shadWeight), 0, 255);
//                int g = glm::clamp(int((colour.green * brightness + s*255.0f) * shadWeight), 0, 255);
//                int b = glm::clamp(int((colour.blue  * brightness + s*255.0f) * shadWeight), 0, 255);
//
//                uint32_t c = (255 << 24) + (r << 16) + (g << 8) + b;
//                window.setPixelColour(u, v, c);
//
//            }
//        }
//    }
//}









//handling the SDL event depending on its type
void handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {

        //moving camera position
       if (event.key.keysym.sym == SDLK_LEFT) {
       std::cout << "LEFT" << std::endl;
       left = true;
       }
       else if (event.key.keysym.sym == SDLK_RIGHT) {
       std::cout << "RIGHT" << std::endl;
       right = true;
       }
       else if (event.key.keysym.sym == SDLK_UP) {
       std::cout << "UP" << std::endl;
       up = true;
       }
       else if (event.key.keysym.sym == SDLK_DOWN) {
       std::cout << "DOWN" << std::endl;
       down = true;
       }
       else if (event.key.keysym.sym == SDLK_q) {
       std::cout << "IN" << std::endl;
       zoomIn = true;
       }
       else if (event.key.keysym.sym == SDLK_e) {
       std::cout << "OUT" << std::endl;
       zoomOut = true;
       }


        //Some kind of texture toggle keypress in the future
       else if (event.key.keysym.sym == SDLK_t) {
            std::cout << "t" << std::endl;
            text = true;
        }


        //orbiting keypresses
        else if (event.key.keysym.sym == SDLK_w) {
            std::cout << "w" << std::endl;
            upOrb = true;
        }
        else if (event.key.keysym.sym == SDLK_a) {
            std::cout << "a" << std::endl;
            leftOrb = true;
        }
        else if (event.key.keysym.sym == SDLK_d) {
            std::cout << "d" << std::endl;
            rightOrb = true;
        }
       else if (event.key.keysym.sym == SDLK_s) {
            std::cout << "s" << std::endl;
            downOrb = true;
        }

        //now switching between different modes
       else if (event.key.keysym.sym == SDLK_b) {
            std::cout << "b" << std::endl;
            wireF = true;
            rast = false;
            ray = false;
        }
       else if (event.key.keysym.sym == SDLK_n) {
            std::cout << "s" << std::endl;
            rast = true;
            ray = false;
            wireF = false;

        }
       else if (event.key.keysym.sym == SDLK_m) {
            std::cout << "m" << std::endl;
            ray = true;
            wireF = false;
            rast = false;

        }

        //automatic orbit conditional on keypress to turn on and off
        else if (event.key.keysym.sym == SDLK_o) {
            if (orbitVar == false) {
                 std::cout << "o" << std::endl;
                 orbitVar = true;
            }
            else {
                std::cout << "o" << std::endl;
                orbitVar = false;
            }
        }

        else if (event.key.keysym.sym == SDLK_ESCAPE) { // Press Esc to exit
            std::cout << "Exiting yeehaw..." << std::endl;
            running = false;
        }

        else if (event.key.keysym.sym == SDLK_r) { // Press Esc to exit
            if (mirror == false) {
                 std::cout << "mirror on" << std::endl;
                 mirror = true;
            }
            else {
                std::cout << "mirror off" << std::endl;
                mirror = false;
            }
        }


    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
       window.savePPM("output.ppm");
       window.saveBMP("output.bmp");
    }
}




//reset the depthBuffer to appropriate - check where this is called to try and solve the occlusion problem properly!
void depthBufReset(std::array<std::array<float, WIDTH>, HEIGHT> &depthBuffer) {
    for (int y = 0; y < HEIGHT; y++)
        for (int x = 0; x < WIDTH; x++)
            depthBuffer[y][x] = 0.0f;
}




//general render function which changes based on the variables used
void render(const std::vector<ModelTriangle> &triVec, DrawingWindow &window, glm::vec3 &camPos, float focalLength, std::array<std::array<float, WIDTH>, HEIGHT> &depthBuffer, glm::vec3 lightSource) {
    if (wireF) {
        window.clearPixels();
        depthBufReset(depthBuffer);
        wireFrameRender(triVec, window, camPos, focalLength, depthBuffer);
    }
    else if (rast) {
        window.clearPixels();
        depthBufReset(depthBuffer);
        rasterRender(triVec, window, camPos, focalLength, depthBuffer);
    }
    else if (ray) {
        window.clearPixels();
        rayTraceRender(triVec, window, camPos, focalLength, lightSource);
    }
}



//main function takes a pointer to argument array?
int main(int argc, char *argv[]) {
    srand(time(0));

    // Initialise important variables - adjustable
    glm::vec3 camPos(0.0f, 0.0f, 4.0f);  // Camera in front of the box
    float focalLength = 320.0f;          // Scale the image
    oriMat = glm::mat3(1.0f);  // Identity orientation matrix to begin with
    glm::vec3 lightSource(0.0f, 0.7f, 1.0f); //position of the source of light
    //now create a light source made up of multiple vec3s
    std::vector<glm::vec3> multiLight;

    std::array<std::array<float, WIDTH>, HEIGHT> depthBuffer; // depth buffer
    depthBufReset(depthBuffer);


    // Load the obj and mtl files for cornell box and sphere
//    triVec = readObjFile("cornell-box.obj", 0.35f, "cornell-box.mtl");
    std::vector<ModelTriangle> sVec = readObjFile("sphere.obj", 0.35f, "sphere.mtl");
    triVec = readObjFile("textured-cornell-box.obj", 0.35f, "textured-cornell-box.mtl");
    triVec.insert(triVec.end(), sVec.begin(), sVec.end());

//    chooseObjects(true, false);

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
        render(triVec, window, camPos, focalLength, depthBuffer, lightSource);

        // TRANSLATION
        if (left) {
            camPos = translatePos(camPos, -move, 0, 0);
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

        //FLY THROUGH + ORBIT? Do a cool camera sequence


        //mirror
         makeTriReflective(triVec, mirror);
         //texture
//         chooseObjects(true, text);


        // Render frame
        window.renderFrame();
    }
}