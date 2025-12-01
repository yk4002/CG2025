#pragma once

#include <string>
#include <vector>
#include <glm/glm.hpp>

std::vector<std::string> split(const std::string &line, char delimiter);
glm::vec3 baryCoords(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 r);
