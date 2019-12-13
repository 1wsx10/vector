#pragma once
#include <iostream>

#ifndef PI
// nasa approved 15 digits of PI
#define PI 3.141592653589793
#endif

#define IS_APPROX_0(X, DELTA) (fabs(X) < (DELTA))

class v3d {
	public:
		static const v3d zero;
		static const v3d unit;
		static const v3d X;
		static const v3d Y;
		static const v3d Z;

		double x;
		double y;
		double z;

		v3d(void);
		//constexpr v3d(double x, double y, double z);
		constexpr v3d(double x, double y, double z): x(x), y(y), z(z) {}

		// static methods here do not modify the inputs
		// non-static versions modify the class its called on

		void set(double x, double y, double z);

		static double angle(const v3d&, const v3d&);

		double dot(const v3d&) const;
		static double dot(const v3d&, const v3d&);

		v3d& cross(const v3d&);
		static v3d cross(v3d, const v3d&);

		double distance(const v3d&) const;
		static double distance(const v3d&, const v3d&);

		double length_squared() const;
		static double length_squared(const v3d&);

		double length() const;
		static double length(const v3d&);

		v3d& normalise();
		static v3d normalise(v3d);

		v3d& try_normalise();
		static v3d try_normalise(v3d);

		v3d &project(const v3d&);
		static v3d project(v3d, const v3d&);

		v3d &reject(const v3d&);
		static v3d reject(v3d, const v3d&);

		v3d &rotate(double angle, const v3d &axis);
		static v3d rotate(v3d to_rotate, double angle, const v3d &axis);

		v3d &mirror(const v3d &axis);
		static v3d mirror(v3d to_mirror, const v3d &axis);

		// * operator is not canonical, so this is explicit how it works
		v3d &mul_elems(const v3d &axis);
		static v3d mul_elems(v3d to_mirror, const v3d &axis);


		bool is_valid() const;
		static bool is_valid(const v3d&);

		bool equals(const v3d&, const double&) const;
		static bool equals(const v3d&, const v3d&, const double&);

		// easy array-like access
		operator const double*() const;
		operator double*();

		v3d operator+(const v3d&) const;
		v3d operator-(const v3d&) const;

		v3d operator-() const;// a = -my_v3d;

		v3d& operator+=(const v3d&);
		v3d& operator-=(const v3d&);

		bool operator==(const v3d&) const;
		bool operator!=(const v3d&) const;

		v3d operator*(double) const;
		v3d operator/(double) const;
		friend v3d operator*(double, const v3d&);
		friend v3d operator/(double, const v3d&);

		v3d& operator*=(double);
		v3d& operator/=(double);

		friend std::ostream& operator<<(std::ostream&, const v3d&);

		static constexpr bool ensure_contiguous_data();
};



constexpr bool v3d::ensure_contiguous_data() {
	v3d val = {1,1,1};

	// check x is at the start
	if((void*)&val.x != (void*)&val)
		return false;

	// check y is after x
	if(&val.x+1 != &val.y)
		return false;

	// check z is after y
	if(&val.y+1 != &val.z)
		return false;

	// check size is 3 doubles
	if(sizeof(v3d) != sizeof(double)*3)
		return false;

	// all tests passed
	return true;
}


