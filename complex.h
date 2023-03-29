#include <cstring>
using namespace std;

template <class T>
class complex {
    T real;
    T imga;
    public:
    complex() { real = 0.; imga = 0.; }                                  // constructors
    complex(T r,T i) { real = r; imga = i; }
    complex(T angle) { real = cos(angle); imga = sin(angle); } 
    void operator=(complex<T> c) {
        real = c.real;
	imga = c.imga;
    }
    complex operator+(complex<T> c) {
        complex<T> cc;
	cc.real = real + c.real;
	cc.imga = imga + c.imga;
	return cc;
    }
    complex operator-(complex<T> c) {
        complex<T> cc;
	cc.real = real - c.real;
	cc.imga = imga - c.imga;
	return cc;
    }
    complex operator*(complex<T> c) {
        complex<T> cc;
	cc.real = real*c.real - imga*c.imga;
	cc.imga = real*c.imga + imga*c.real;
	return cc;
    }
    complex operator*(T c) {
        complex<T> cc;
	cc.real = real*c;
	cc.imga = imga*c;
	return cc;
    }
    complex operator/(complex<T> c) {
        complex<T> cc;
	T a = c.real*c.real+c.imga*c.imga;
	cc.real = (real*c.real+imga*c.imga)/a;
	cc.imga = (imga*c.real-real*c.imga)/a;
	return cc;
    }
    complex conjugate() {
        complex<T> cc;
        cc.real = real;
	cc.imga = -1.*imga;
	return cc;
    }
    void operator+=(complex<T> c) {
	real += c.real;
	imga += c.imga;
    }
    void operator-=(complex<T> c) {
	real -= c.real;
	imga -= c.imga;
    }
    void operator*=(complex<T> c) {
        T a,b;
	a = real*c.real - imga*c.imga;
	b = real*c.imga + imga*c.real;
	real = a;
	imga = b;
    }
    void operator*=(T c) {
        real *= c;
	imga *= c;
    }
    void operator/=(complex<T> c) {
        T a,u,v;
	a = c.real*c.real+c.imga*c.imga;
        u = (real*c.real+imga*c.imga)/a;
	v = (imga*c.real-real*c.imga)/a;
	real = u;
	imga = v;
    }
    void setzero() {
        real = 0.;
	imga = 0.;
    }
    void setrealimga(T r,T i) {
        real = r;
	imga = i;
    }
    void setreal(T r) { real = r; }
    void setimga(T i) { imga = i; }
    void setangle(T angle) { real = cos(angle); imga = sin(angle); }
    T realpart() {
        return real;
    }
    T imgapart() {
        return imga;
    }
    void print() {
        cout << real << " " << imga << endl;
    }
};

