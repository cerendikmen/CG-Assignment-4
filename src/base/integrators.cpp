#include "integrators.hpp"

#include "utility.hpp"
#include "particle_systems.hpp"

void eulerStep(ParticleSystem& ps, float step) {
	// YOUR CODE HERE (R1)
	// Implement an Euler integrator.
	State initial;
	State final;
	State f;

	initial = ps.state();
	f = ps.evalF(initial);
	for (int i = 0; i < initial.size(); i++)
	{
		final.push_back(initial[i] + step * f[i]);

	}
	ps.set_state(final);


};

void trapezoidStep(ParticleSystem& ps, float step) {
	// YOUR CODE HERE (R3)
	// Implement a trapezoid integrator.

	State initial;
	State final;
	State final2;
	State f;
	State f2;

	initial = ps.state();
	f = ps.evalF(initial);
	for (int i = 0; i < initial.size(); i++)
	{
		final.push_back(initial[i] + step * f[i]);

	}

	ps.set_state(final);
	f2 = ps.evalF(ps.state());
	for (int i = 0; i < initial.size(); i++)
	{
		final2.push_back(initial[i] + (step/2) * (f[i] + f2[i]));

	}

	ps.set_state(final2);
}

void midpointStep(ParticleSystem& ps, float step) {
	const auto& x0 = ps.state();
	auto n = x0.size();
	auto f0 = ps.evalF(x0);
	auto xm = State(n), x1 = State(n);
	for (auto i = 0u; i < n; ++i) {
		xm[i] = x0[i] + (0.5f * step) * f0[i];
	}
	auto fm = ps.evalF(xm);
	for (auto i = 0u; i < n; ++i) {
		x1[i] = x0[i] + step * fm[i];
	}
	ps.set_state(x1);
}

void rk4Step(ParticleSystem& ps, float step) {
	// EXTRA: Implement the RK4 Runge-Kutta integrator.
}
 
