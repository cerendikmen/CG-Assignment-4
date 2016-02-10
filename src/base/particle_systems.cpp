#include "particle_systems.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>

using namespace std;
using namespace FW;

namespace {

inline Vec3f fGravity(float mass) {
	return Vec3f(0, -9.8f * mass, 0);
}

// force acting on particle at pos1 due to spring attached to pos2 at the other end
inline Vec3f fSpring(const Vec3f& pos1, const Vec3f& pos2, float k, float rest_length) {
	// YOUR CODE HERE (R2)
	// Note that you can visualize the forces using the draw_lines function found in utility.hpp. This might make debugging easier for you.
	float distance = FW::sqrt(FW::pow((pos1.x - pos2.x), 2) + FW::pow((pos1.y - pos2.y), 2) + FW::pow((pos1.z - pos2.z), 2));
	Vec3f spring_force = k*(rest_length - distance)*((pos1 - pos2) / distance);
	return spring_force;
}

inline Vec3f fDrag(const Vec3f& v, float k) {
	// YOUR CODE HERE (R2)
	Vec3f drag_force;
	drag_force = (-1)*k*v;
	return drag_force;
}

} // namespace

void SimpleSystem::reset() {
	state_ = State(1, Vec3f(0, radius_, 0));
}

State SimpleSystem::evalF(const State& state) const {
	State f(1, Vec3f(-state[0].y, state[0].x, 0));
	return f;
}

Points SimpleSystem::getPoints() {
	return Points(1, state_[0]);
}

Lines SimpleSystem::getLines() {
	static const auto n_lines = 50u;
	auto l = Lines(n_lines * 2);
	const auto angle_incr = 2*FW_PI/n_lines;
	for (auto i = 0u; i < n_lines; ++i) {
		l[2*i] = l[2*i+1] =
			Vec3f(radius_ * FW::sin(angle_incr * i), radius_ * FW::cos(angle_incr * i), 0);
	}
	rotate(l.begin(), l.begin()+1, l.end());
	return l;
}

void SpringSystem::reset() {
	const auto start_pos = Vec3f(0.1f, -0.5f, 0.0f);
	const auto spring_k = 30.0f;
	state_ = State(4);
	// YOUR CODE HERE (R2)
	// Set the initial state for a particle system with one particle fixed
	// at origin and another particle hanging off the first one with a spring.
	// Place the second particle initially at start_pos.
	const auto origin_pos = Vec3f(0.0f, 0.0f, 0.0f);
	auto rest_length = FW::sqrt(FW::pow((start_pos.x - origin_pos.x), 2) + FW::pow((start_pos.y - origin_pos.y), 2) + FW::pow((start_pos.z - origin_pos.z), 2));
	spring_ = Spring(0, 1, spring_k, rest_length);
	state_[0] = origin_pos;
	state_[1] = Vec3f(0.0f, 0.0f, 0.0f);
	state_[2] = start_pos;
	state_[3] = Vec3f(0.0f, 0.0f, 0.0f);

}

State SpringSystem::evalF(const State& state) const {
	const auto drag_k = 0.5f;
	const auto mass = 1.0f;
	State f(4);
	// YOUR CODE HERE (R2)
	// Return a derivative for the system as if it was in state "state".
	// You can use the fGravity, fDrag and fSpring helper functions for the forces.
	Vec3f dragF;
	Vec3f springF;
	Vec3f gravityF;
	Vec3f force;
	f[0] = state[1];
	f[1] = 0;
	dragF = fDrag(state[3], drag_k);
	springF = fSpring(state[2], state[0], spring_.k, spring_.rlen);
	gravityF = fGravity(mass);
	force = dragF + springF + gravityF;
	f[2] = state[3];
	f[3] = force;
	return f;
}

Points SpringSystem::getPoints() {
	auto p = Points(2);
	p[0] = state_[0]; p[1] = state_[2];
	return p;
}

Lines SpringSystem::getLines() {
	auto l = Lines(2);
	l[0] = state_[0]; l[1] = state_[2];
	return l;
}


void PendulumSystem::reset() {
	const auto spring_k = 1000.0f;
	const auto start_point = Vec3f(0);
	const auto end_point = Vec3f(0.05, -1.5, 0);
	state_ = State(2*n_);
	int i;
	// YOUR CODE HERE (R4)
	// Set the initial state for a pendulum system with n_ particles
	// connected with springs into a chain from start_point to end_point.
	float scale_x = (end_point.x - start_point.x) / (n_ - 1);
	float scale_y = (end_point.y - start_point.y) / (n_ - 1);
	float rest_length = FW::sqrt(FW::pow((end_point.x - start_point.x), 2) + FW::pow((end_point.y - start_point.y), 2) + FW::pow((end_point.z - start_point.z), 2));
	float piece_length = rest_length / (n_ - 1);
	for ( i = 0; i < n_-1; i++)
	{
		state_[2*i] = Vec3f(start_point.x + i * scale_x , start_point.y + i * scale_y, 0);
		state_[2*i + 1] = Vec3f(0.0f, 0.0f, 0.0f);
		springs_.push_back( Spring(i, i + 1, spring_k,piece_length));
	}
	state_[2 * i] = end_point;
	state_[2 * i + 1] = Vec3f(0.0f, 0.0f, 0.0f);
}
  
State PendulumSystem::evalF(const State& state) const {
	const auto drag_k = 0.5f;
	const auto mass = 0.5f;
	auto f = State(2*n_);
	// YOUR CODE HERE (R4)
	// As in R2, return a derivative of the system state "state".
	Vec3f dragF;
	Vec3f first_springF;
	Vec3f second_springF;
	Vec3f gravityF;
	Vec3f force;
	gravityF = fGravity(mass);
	f[0] = state[1];
	f[1] = 0;

	for (int i = 1; i < n_ -1; i++)
	{
		f[2 * i] = state[2 * i + 1];
		dragF = fDrag(state[2*i + 1], drag_k);
		first_springF = fSpring(state[2*i], state[2*i-2], springs_[i-1].k, springs_[i-1].rlen);
		second_springF = fSpring(state[2 * i], state[2 * i + 2], springs_[i].k, springs_[i].rlen);
		force = dragF + gravityF + first_springF + second_springF;
		f[2 * i + 1] = force;
	}
	f[2 * (n_ - 1)] = state[2 * (n_ - 1) + 1];
	dragF = fDrag(state[2 * (n_ - 1) + 1], drag_k);
	first_springF = fSpring(state[2 * (n_ - 1)], state[2 * (n_ - 1) - 2], springs_[(n_ - 1) - 1].k, springs_[(n_ - 1) - 1].rlen);
	force = dragF + gravityF + first_springF;
	f[2 * (n_ - 1) + 1] = force;

	return f;
}

Points PendulumSystem::getPoints() {
	auto p = Points(n_);
	for (auto i = 0u; i < n_; ++i) {
		p[i] = state_[i*2];
	}
	return p;
}

Lines PendulumSystem::getLines() {
	auto l = Lines();
	for (const auto& s : springs_) {
		l.push_back(state_[2*s.i1]);
		l.push_back(state_[2*s.i2]);
	}
	return l;
}

void ClothSystem::reset() {
	const auto spring_k = 300.0f;
	const auto width = 1.5f, height = 1.5f; // width and height of the whole grid
	state_ = State(2 * x_*y_);
	// YOUR CODE HERE (R5)
	// Construct a particle system with a x_ * y_ grid of particles,
	// connected with a variety of springs as described in the handout:
	// structural springs, shear springs and flex springs.
	for (int i = 0; i < x_; i++){
		for (int j = 0; j < y_; j++)
		{
			float x = (float)i / ((float)x_ - 1.0f) * height - height / 2.0f;
			float y = (float)j / ((float)y_ - 1.0f) * width * -1.0f;
			int pos = 2 * i*x_ + 2 * j;
			state_[pos] = Vec3f(x, 0, y);
		}
	}
	if (springs_.size() == 0){
		float restl = 1.5f / 19.0f;
		for (int i = 0; i < x_ * y_; i++){
			if ((i + 1) % 20 != 0){
				Spring k(i, i + 1, spring_k, restl);
				springs_.push_back(k);
			}
			if (i + 20 < 400){
				Spring l(i, i + 20, spring_k, restl);
				springs_.push_back(l);
			}
			if ((i - 20) > 0 && (i + 1) % 20 != 0){
				Spring m(i, i - 19, spring_k, restl * sqrt(2));
				springs_.push_back(m);
			}
			if (i + 20 < 400 && (i + 1) % 20 != 0){
				Spring m(i, i + 21, spring_k, restl * sqrt(2));
				springs_.push_back(m);
			}
		}
	}
}
State ClothSystem::evalF(const State& state) const {
	const auto drag_k = 0.08f;
	const auto n = x_ * y_;
	static const auto mass = 0.025f;
	auto f = State(2 * n);
	// YOUR CODE HERE (R5)
	// This will be much like in R2 and R4.
	//float cross = state[0].x * state[0].x
	for (int i = 1; i < n; i++){
		if (i == 380)
			continue;
		f[2 * i] = state[2 * i + 1];
		Vec3f kv = fGravity(mass);
		Vec3f drg = fDrag(state[2 * i + 1], drag_k);
		Vec3f fspr(0);
		
		if ((i + 1) % 20 != 0)
			fspr += fSpring(state[2 * i], state[2 * (i + 1)], springs_[0].k, springs_[0].rlen);
		if ((i + 20) < 400)
			fspr += fSpring(state[2 * i], state[2 * (i + 20)], springs_[0].k, springs_[0].rlen);
		if (i % 20 != 0)
			fspr += fSpring(state[2 * i], state[2 * (i - 1)], springs_[0].k, springs_[0].rlen);
		if ((i - 20) >= 0)
			fspr += fSpring(state[2 * i], state[2 * (i - 20)], springs_[0].k, springs_[0].rlen);
		if ((i + 20) < 400 && (i + 1) % 20 != 0)
			fspr += fSpring(state[2 * i], state[2 * (i + 21)], springs_[0].k, (springs_[0].rlen)* sqrt(2));
		if ((i - 20) > 0 && (i % 20 != 0))
			fspr += fSpring(state[2 * i], state[2 * (i - 21)], springs_[0].k, springs_[0].rlen * sqrt(2));
		if ((i + 20) < 400 && i % 20 != 0)
			fspr += fSpring(state[2 * i], state[2 * (i + 19)], springs_[0].k, springs_[0].rlen * sqrt(2));
		if ((i - 20) > 0 && ((i + 1) % 20 != 0))
			fspr += fSpring(state[2 * i], state[2 * (i - 19)], springs_[0].k, springs_[0].rlen * sqrt(2));

		if ((i + 1) % 20 != 0 && (i + 2) % 20 != 0)
			fspr += fSpring(state[2 * i], state[2 * (i + 2)], springs_[0].k, 2 * springs_[0].rlen);
		if ((i + 40) < 400)
			fspr += fSpring(state[2 * i], state[2 * (i + 40)], springs_[0].k, 2 * springs_[0].rlen);
		if ((i % 20 != 0) && (i % 20 != 1))
			fspr += fSpring(state[2 * i], state[2 * (i - 2)], springs_[0].k, 2 * springs_[0].rlen);
		if ((i - 40) >= 0)
			fspr += fSpring(state[2 * i], state[2 * (i - 40)], springs_[0].k, 2 * springs_[0].rlen);
		f[2 * i + 1] = (kv + drg + fspr) / mass;
	}
	return f;
}

Points ClothSystem::getPoints() {
	auto n = x_ * y_;
	auto p = Points(n);
	for (auto i = 0u; i < n; ++i) {
		p[i] = state_[2*i];
	}
	return p;
}

Lines ClothSystem::getLines() {
	auto l = Lines();
	for (const auto& s : springs_) {
		l.push_back(state_[2*s.i1]);
		l.push_back(state_[2*s.i2]);
	}
	return l;
}
