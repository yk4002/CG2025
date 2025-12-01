#pragma once

#include <glm/glm.hpp>
#include <string>
#include <array>
#include "Colour.h"
#include "TexturePoint.h"
#include <vector>

struct ModelTriangle {
	std::array<glm::vec3, 3> vertices{};
	std::vector<TexturePoint> texturePoints{};
	Colour colour{};
	glm::vec3 normal{};
	bool isMirror = false;

	ModelTriangle();
	ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Colour trigColour);
	friend std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle);
};
