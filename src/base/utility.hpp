#pragma once

#include "../framework/base/Math.hpp"

#include <cassert>
#include <iostream>
#include <vector>

template <typename T>
void printVec3(const T& t) {
	cout << "(" << t.x << " " << t.y << " " << t.z << ")" << endl;
}


// Inline functions to help with drawing
inline void glVertex(const FW::Vec3f& a)
{
	glVertex3fv(a.getPtr());
}

inline void glColor(const FW::Vec3f& a)
{
	glColor3fv(a.getPtr());
}

inline void glNormal(const FW::Vec3f& a)
{
	glNormal3fv(a.getPtr());
}

// This function draws lines and/or points at the given positions.
// If lines is true, line segmens are drawn between adjacent vertices
// in the positions vector (the first point and the second, the third
// and the fourth etc), so if you want a segmented line path A-B-C,
// the vector should be (A,B,B,C). Points of the given size (in pixels)
// will be drawn at the positions. The colors vector should contain a
// color value for each vertex in the positions vector.

// Note that the positions will be relative to the current opengl transformation state.
// In practice this means that this should be called only in App::render after the
// call to camera_.ApplyModelview();.

inline void draw_lines(std::vector<FW::Vec3f> positions, bool lines = true, float point_size = .0f, std::vector<FW::Vec3f> colors = std::vector<FW::Vec3f>()) {

	glPushAttrib(GL_ENABLE_BIT);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glDisable(GL_LIGHTING);

	bool colored = (positions.size() == colors.size());

	glColor(FW::Vec3f(1.0f));

	if (lines) {
		glBegin(GL_LINES);
		for (int i = 0; i<positions.size(); ++i) {
			if (colored)
				glColor(colors[i]);
			glVertex(positions[i]);
		}
		glEnd();
	}

	if (point_size > 0) {
		glPointSize(point_size);
		glBegin(GL_POINTS);
		for (int i = 0; i<positions.size(); ++i) {
			if (colored)
				glColor(colors[i]);
			glVertex(positions[i]);
		}
		glEnd();
	}

	glPopAttrib();
}

inline void glLoadMatrix(const FW::Mat4f& m)
{
	glLoadMatrixf(m.getPtr());
}

inline void glMultMatrix(const FW::Mat4f& m)
{
	glMultMatrixf(m.getPtr());
}

// Axis must be unit length. Angle in radians.
inline FW::Mat4f rotation4f(const FW::Vec3f& axis, float angle)
{
	FW::Mat3f R = FW::Mat3f::rotation(axis, angle);
	// Copy to homogeneous matrix
	FW::Mat4f H;
	H.setIdentity();
	for (unsigned i = 0; i < 3; ++i)
		for (unsigned j = 0; j < 3; ++j)
			H(i, j) = R(i, j);
	return H;
}

inline void printTranspose(const FW::Vec3f& a)
{
	std::cout << a.x << " " << a.y << " " << a.z;
}

inline void printTranspose(const FW::Vec4f& a)
{
	std::cout << a.x << " " << a.y << " " << a.z << " " << a.w;
}


// Construct a matrix from the columns [a, b, c, d].
inline FW::Mat4f makeMat4f(const FW::Vec4f& a, const FW::Vec4f& b, const FW::Vec4f& c, const FW::Vec4f& d) {
	FW::Mat4f A;
	A.col(0) = a;
	A.col(1) = b;
	A.col(2) = c;
	A.col(3) = d;
	return A;
}

// http://stackoverflow.com/questions/2270726/how-to-determine-the-size-of-an-array-of-strings-in-c
template <typename T, std::size_t N>
char (&static_sizeof_array( T(&)[N] ))[N];   // declared, not defined
#define SIZEOF_ARRAY( x ) sizeof(static_sizeof_array(x))

