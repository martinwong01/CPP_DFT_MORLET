#include <cstring>
using namespace std;

template <class Type>
class complex {
    Type real;
    Type imga;
    public:
    complex<Type>() { real = 0.; imga = 0.; }                                  // constructors
    complex<Type>(Type r,Type i) { real = r; imga = i; }
    complex<Type>(Type angle) { real = cos(angle); imga = sin(angle); } 
    void operator=(complex<Type> c) {
        real = c.real;
	imga = c.imga;
    }
    complex<Type> operator+(complex<Type> c) {
        complex<Type> cc;
	cc.real = real + c.real;
	cc.imga = imga + c.imga;
	return cc;
    }
    complex<Type> operator-(complex<Type> c) {
        complex<Type> cc;
	cc.real = real - c.real;
	cc.imga = imga - c.imga;
	return cc;
    }
    complex<Type> operator*(complex<Type> c) {
        complex<Type> cc;
	cc.real = real*c.real - imga*c.imga;
	cc.imga = real*c.imga + imga*c.real;
	return cc;
    }
    complex<Type> operator*(Type c) {
        complex<Type> cc;
	cc.real = real*c;
	cc.imga = imga*c;
	return cc;
    }
    complex<Type> operator/(complex<Type> c) {
        complex<Type> cc;
	Type a = c.real*c.real+c.imga*c.imga;
	cc.real = (real*c.real+imga*c.imga)/a;
	cc.imga = (imga*c.real-real*c.imga)/a;
	return cc;
    }
    complex<Type> turnleft() {
        complex<Type> cc;
	cc.real = -imga;
	cc.imga = real;
	return cc;
    }
    complex<Type> turnright() {
        complex<Type> cc;
	cc.real = imga;
	cc.imga = -real;
	return cc;
    }
    complex<Type> reverse() {
        complex<Type> cc;
	cc.real = -real;
	cc.imga = -imga;
	return cc;
    }
    complex<Type> conjugate() {
        complex<Type> cc;
        cc.real = real;
	cc.imga = -1.*imga;
	return cc;
    }
    void operator+=(complex<Type> c) {
	real += c.real;
	imga += c.imga;
    }
    void operator-=(complex<Type> c) {
	real -= c.real;
	imga -= c.imga;
    }
    void operator*=(complex<Type> c) {
        Type a,b;
	a = real*c.real - imga*c.imga;
	b = real*c.imga + imga*c.real;
	real = a;
	imga = b;
    }
    void operator*=(Type c) {
        real *= c;
	imga *= c;
    }
    void operator/=(complex<Type> c) {
        Type a,u,v;
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
    void setrealimga(Type r,Type i) {
        real = r;
	imga = i;
    }
    void setreal(Type r) { real = r; }
    void setimga(Type i) { imga = i; }
    void setangle(Type angle) { real = cos(angle); imga = sin(angle); }
    Type realpart() {
        return real;
    }
    Type imgapart() {
        return imga;
    }
    void print() {
        cout << real << " " << imga << endl;
    }
};

