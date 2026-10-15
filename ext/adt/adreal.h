// -*- Mode: c++ -*-
#ifndef __ADREAL_H__
#define __ADREAL_H__

#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include <float.h>
#ifdef COUNTER_ON
#include <common/Counter.h>
#endif


// Defines automtic differentiation scalar 
// with forward computation of derivatives of order 1 or order 1 and 2
// supports arithmetic operators +,-,*,/, unary minus, comparison operators ==, !=, <, >, <=, >=
// assignment operators =, +=, -=, *-, /=
// for comparison operators one argument may be of corresponding constant
// scalar type
// for assignment operators rhs can  be of corresponding constant type
// ANSI C math library support: 
// acos(),asin(),atan(), atan2(), cos(), exp(), log(), pow(), sin(), sqrt(), tan(), fabs()
// in addition defines max and min
// TODO: log10(), sinh(), cosh(), tanh() 
//
// The following standard math functions  are not supported
// as it is not clear if anything meaningful can or needs to be done 
// ceil(), floor(), fmod(), frexp(), ldexp(), modf()


typedef unsigned int uint;

using namespace std;

// the "right" way to do this is to make expressions into functors etc 
// but it is not clear how good is the optimizer is at unwinding it all
// so using a macro
// iterate over all gradient and hessian entries 
// assigning the function, gradient and hessian  values of expressions passed as parameters


#define opinst(val_expr, grad_expr, hess_expr) \
  adreal<NUM_VARS,DO_HESS,constreal> temp; \
  temp.value() = (val_expr); \
  for(int i = 0; i < NUM_VARS; i++) { \
    temp.gradient(i) = (grad_expr); \
    for(int j = i; j < NUM_VARS*DO_HESS; j++) {\
      temp.hessian(i,j) = (hess_expr); \
      temp.hessian(j,i) = temp.hessian(i,j); \
        } \
  } \
  return temp

#define sqr(x) ((x)*(x))
#define cub(x) ((x)*(x)*(x))

#include "fixedarray.h"

// the main AD scalar class; for each variable, there is an array
// of NUM_VARS entries of the gradient and NUM_VARS^2 hessian matrix
// DO_HESS has to be 0 or 1, otherwise all memory management breaks; 
// all constructors must check this

template <int NUM_VARS, int DO_HESS, class constreal>
class adreal {

public: 

    adreal() {     assert(DO_HESS == 0 || DO_HESS == 1); }


    adreal(constreal v, const constreal g[NUM_VARS]=0,  const constreal h[NUM_VARS]=0) { 
        assert(DO_HESS == 0 || DO_HESS == 1);
        val = v;
        if(g) memcpy(derivdata, g, sizeof(constreal)*NUM_VARS); else memset(derivdata, 0, sizeof(constreal)*NUM_VARS); 
        if(h) memcpy(&(derivdata[NUM_VARS]),h, sizeof(constreal)*DO_HESS*NUM_VARS*NUM_VARS); else memset(&(derivdata[NUM_VARS]),0, sizeof(constreal)*DO_HESS*NUM_VARS*NUM_VARS); 
    }

    ~adreal() {}
  
    // iniitalizes the variable as an independent variable with index i and
    // value v
    // the gradient and hessian for the dependent variables 
    // will have derivatives with respect to this variable in positions i 
    // and row/column i respectively

    void set_independent( constreal v, uint n) { 
        assert(n < NUM_VARS);
        val = v;    
        memset(derivdata, 0, sizeof(constreal)*NUM_VARS);  
        memset(&(derivdata[NUM_VARS]), 0, sizeof(constreal)*DO_HESS*NUM_VARS*NUM_VARS);
        gradient(n) = 1; 
    }


    // assignment and copy constructors
    adreal& operator= (const adreal& a) {  memcpy(this,&a, sizeof(adreal<NUM_VARS,DO_HESS,constreal>));  return *this;  }


    adreal(const adreal& a) { memcpy(this, &a, sizeof(adreal<NUM_VARS,DO_HESS,constreal>)); }

    // accessors
    constreal    value()    const { return val; }
    constreal&   value()   { return val; }

    constreal    gradient(uint i) const { assert(i < NUM_VARS); return derivdata[i]; } 
    constreal&   gradient(uint i)       { assert(i < NUM_VARS); return derivdata[i]; } 
    GradientType<NUM_VARS,constreal>  gradient() const {  return *( (GradientType<NUM_VARS,constreal>*)(derivdata));   }

    constreal    hessian(uint i, uint j) const { assert(i < NUM_VARS && j < NUM_VARS); return derivdata[NUM_VARS+i*NUM_VARS+j]; } 
    constreal&   hessian(uint i, uint j)       { assert(i < NUM_VARS && j < NUM_VARS); return derivdata[NUM_VARS+i*NUM_VARS+j]; } 
    HessianType<NUM_VARS,constreal> hessian()  const { assert(DO_HESS); return *( (HessianType<NUM_VARS,constreal>*)(&derivdata[NUM_VARS])); }

  
    // assignments 
    // certainly can be optimized, but adding them in simplest form for
    // completeness now 
    adreal& operator+=(const adreal& ad  ) { return (*this) = (*this)+ad; }
    adreal& operator-=(const adreal& ad  ) { return (*this) = (*this)-ad; }
    adreal& operator*=(const adreal& ad  ) { return (*this) = (*this)*ad; }
    adreal& operator/=(const adreal& ad  ) { return (*this) = (*this)/ad; }
    adreal& operator+=(const constreal& c) { return (*this) = (*this)+c;  }
    adreal& operator-=(const constreal& c) { return (*this) = (*this)-c;  }
    adreal& operator*=(const constreal& c) { return (*this) = (*this)*c;  }
    adreal& operator/=(const constreal& c) { return (*this) = (*this)/c;  }

    // exceptions: the only operators defined  by hand: operator/ except 1/c
    // for some reason hessian expressions are a pain to generate in Maple in
    // reasonably optimal form
    // so for now just reduce it to multiplication
    // extra temp overhead
    adreal operator/(const adreal& a) const { return (*this)*(1.0/a); }
    adreal operator/ (constreal c) const { constreal f = 1.0/c; return f*(*this); }
    // unary minus
    adreal operator-(void) const { opinst(-val, -gradient(i), -hessian(i,j)); }
    
  
private:
    constreal val;
    // all derivative first (and possibly second) data goes in here
    constreal derivdata[NUM_VARS + DO_HESS*NUM_VARS*NUM_VARS];
};

template <int NUM_VARS,int DO_HESS, class constreal>
inline bool operator< (const adreal<NUM_VARS,DO_HESS,constreal>& ad1, const adreal<NUM_VARS,DO_HESS,constreal>& ad2) { return (ad1.value() < ad2.value()); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator<=(const adreal<NUM_VARS,DO_HESS,constreal>& ad1, const adreal<NUM_VARS,DO_HESS,constreal>& ad2) { return (ad1.value() <=ad2.value()); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator> (const adreal<NUM_VARS,DO_HESS,constreal>& ad1, const adreal<NUM_VARS,DO_HESS,constreal>& ad2) { return (ad1.value() > ad2.value()); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator>=(const adreal<NUM_VARS,DO_HESS,constreal>& ad1, const adreal<NUM_VARS,DO_HESS,constreal>& ad2) { return (ad1.value() >=ad2.value()); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator==(const adreal<NUM_VARS,DO_HESS,constreal>& ad1, const adreal<NUM_VARS,DO_HESS,constreal>& ad2) { return (ad1.value() ==ad2.value()); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator!=(const adreal<NUM_VARS,DO_HESS,constreal>& ad1, const adreal<NUM_VARS,DO_HESS,constreal>& ad2) { return (ad1.value() !=ad2.value()); }

template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator< (const constreal c, const adreal<NUM_VARS,DO_HESS,constreal>& ad)  { return (c < ad.value()); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator<=(const constreal c, const adreal<NUM_VARS,DO_HESS,constreal>& ad)  { return (c <=ad.value()); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator> (const constreal c, const adreal<NUM_VARS,DO_HESS,constreal>& ad)  { return (c > ad.value()); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator>=(const constreal c, const adreal<NUM_VARS,DO_HESS,constreal>& ad)  { return (c >=ad.value()); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator==(const constreal c, const adreal<NUM_VARS,DO_HESS,constreal>& ad)  { return (c ==ad.value()); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator!=(const constreal c, const adreal<NUM_VARS,DO_HESS,constreal>& ad)  { return (c !=ad.value()); }

template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator< (const adreal<NUM_VARS,DO_HESS,constreal>& ad, const constreal c)  { return (ad.value() < c); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator<=(const adreal<NUM_VARS,DO_HESS,constreal>& ad, const constreal c)  { return (ad.value() <=c); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator> (const adreal<NUM_VARS,DO_HESS,constreal>& ad, const constreal c)  { return (ad.value() > c); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator>=(const adreal<NUM_VARS,DO_HESS,constreal>& ad, const constreal c)  { return (ad.value() >=c); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator==(const adreal<NUM_VARS,DO_HESS,constreal>& ad, const constreal c)  { return (ad.value() ==c); }
template <int NUM_VARS,int DO_HESS, class constreal> 
inline bool operator!=(const adreal<NUM_VARS,DO_HESS,constreal>& ad, const constreal c)  { return (ad.value() !=c); }


// functions with derivative discontinuities
// derivative at the discontinuity is defined by arbitrary choice of 
// one-sided limit
// defining min/max consistently in terms of fabs

template <int NUM_VARS,int DO_HESS, class constreal> 
adreal<NUM_VARS,DO_HESS,constreal>  fabs(const adreal<NUM_VARS,DO_HESS,constreal>& ad) {
    return ((ad >= constreal(0) )? ad :(-ad));
}

template <int NUM_VARS,int DO_HESS, class constreal> 
adreal<NUM_VARS,DO_HESS,constreal>  max(const adreal<NUM_VARS,DO_HESS,constreal>& ad1,const adreal<NUM_VARS,DO_HESS,constreal>& ad2) {
    return 0.5*(fabs(ad1+ad2) + fabs(ad1-ad2));
}

template <int NUM_VARS,int DO_HESS, class constreal> 
adreal<NUM_VARS,DO_HESS,constreal>  max(const adreal<NUM_VARS,DO_HESS,constreal>& ad1,const constreal& c2) {
    return 0.5*(fabs(ad1+c2) + fabs(ad1-c2));
}

template <int NUM_VARS,int DO_HESS, class constreal> 
adreal<NUM_VARS,DO_HESS,constreal>  max(const constreal& c1,const adreal<NUM_VARS,DO_HESS,constreal>& ad2) {
    return 0.5*(fabs(c1+ad2) + fabs(c1-ad2));
}

template <int NUM_VARS,int DO_HESS, class constreal> 
adreal<NUM_VARS,DO_HESS,constreal>  min(const adreal<NUM_VARS,DO_HESS,constreal>& ad1,const adreal<NUM_VARS,DO_HESS,constreal>& ad2) {
    return 0.5*(fabs(ad1+ad2) - fabs(ad1-ad2));
}

template <int NUM_VARS,int DO_HESS, class constreal> 
adreal<NUM_VARS,DO_HESS,constreal>  min(const adreal<NUM_VARS,DO_HESS,constreal>& ad1,const constreal& c2) {
    return 0.5*(fabs(ad1+c2) - fabs(ad1-c2));
}

template <int NUM_VARS,int DO_HESS, class constreal> 
adreal<NUM_VARS,DO_HESS,constreal>  min(const constreal& c1,const adreal<NUM_VARS,DO_HESS,constreal>& ad2) {
    return 0.5*(fabs(c1+ad2) - fabs(c1-ad2));
}



template <int NUM_VARS,int DO_HESS, class constreal> 
ostream& operator<< (ostream& os, const adreal<NUM_VARS,DO_HESS,constreal>& a) {
    os << a.value() <<  " " << a.gradient();
    if(DO_HESS) 
        os << " " << a.hessian();
    os << endl;
    return os;

}

// not defining this as semantics is unclear
//template <int NUM_VARS,int DO_HESS, class constreal> 
//  ostream& operator<< (ostream& os, const adreal<NUM_VARS,DO_HESS,constreal>& a)



// this is used by automatically generated tests
template <int NUM_VARS,int DO_HESS, class constreal> void compareAD
(const adreal<NUM_VARS,DO_HESS,constreal>& ad1, const adreal<NUM_VARS,DO_HESS,constreal>& ad2) { 
    constreal diffv,normv, diffg,normg, diffh,normh;
    normv = fabs(ad1.value());
    diffv = fabs(ad1.value()-ad2.value());

    if( diffv/normv > 1 ) cout << "WEIRD: " << setprecision(20) << 
                              ad1.value() << " " << ad2.value() << " " << diffv << 
                              " " << normv  << endl;
    if(normv > DBL_EPSILON) diffv /= normv;

    normg = 0; diffg = 0;
    for(int i =0; i < NUM_VARS; i++) {
        normg += sqr(ad1.gradient(i));
        diffg += sqr(ad1.gradient(i)-ad2.gradient(i));
    }
    normg = sqrt(normg);
    diffg = sqrt(diffg);
    if(normg > DBL_EPSILON) diffg /= normg;

    if( DO_HESS) { 
        normh = 0; diffh = 0;
        for(int i =0; i < NUM_VARS; i++) {
            for(int j =0; j < NUM_VARS; j++) {
                normh += sqr(ad1.hessian(i,j));
                diffh += sqr(ad1.hessian(i,j)-ad2.hessian(i,j));
            }
        }
        normh = sqrt(normh);
        diffh = sqrt(diffh);
        if(normh > DBL_EPSILON) diffh /= normh;
    }
    else diffh = 0; 
}

// synonym for the member function of the same name
// convenient to have this so that one can have 
// templates that work for ADReal and double the same way

template <int NUM_VARS,int DO_HESS, class constreal>
void set_independent(adreal<NUM_VARS,DO_HESS,constreal>& a, constreal v, uint n) { 
  a.set_independent(v,n);
}

// dummies for float and double 

inline void set_independent(double& a, double v, uint n) { 
  a = v;
}

inline void set_independent(float& a, float v, uint n) { 
  a = v;
}


// automatically generated code include
#include "adreal_maple.h"

#endif
