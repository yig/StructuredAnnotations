#ifndef __FIXED_ARRAY_H__
#define __FIXED_ARRAY_H__

#include <iostream>
using namespace std;
// structures to pass  Gradient and Hessian
template <int NUM_VARS, class constreal> 
class HessianType { 
public:
  static const int NUM_DOF = NUM_VARS;
  constreal& operator()(unsigned int i, unsigned int j) { assert( i < NUM_VARS && j < NUM_VARS); return hess[i][j]; }
  const constreal& operator()(unsigned int i, unsigned int j) const {    assert( i < NUM_VARS && j < NUM_VARS); return hess[i][j]; }  


private:
  constreal hess[NUM_VARS][NUM_VARS];
};

template <int NUM_VARS, class constreal> 
class GradientType { 
public:
  static const int NUM_DOF = NUM_VARS;
  constreal& operator()(unsigned int i) { assert( i < NUM_VARS);  return grad[i]; }  
  const constreal& operator()(unsigned int i) const { assert( i < NUM_VARS); return grad[i]; }

  constreal& operator[](unsigned int i) { assert( i < NUM_VARS);  return grad[i]; }  
  const constreal& operator[](unsigned int i) const { assert( i < NUM_VARS); return grad[i]; }
  constreal* data() { return grad; }
  const constreal* data() const { return grad; }

private:
  constreal grad[NUM_VARS];    
};


template <int NUM_VARS, class constreal> 
ostream& operator<<(ostream& os, const GradientType<NUM_VARS,constreal>& g) {
    os << "[";
    for(int j = 0; j < NUM_VARS; j++)   os << g(j) << ( (j == NUM_VARS-1)?"]\n":", "); 
    return os;
}

template <int NUM_VARS, class constreal>
ostream& operator<<(ostream& os, const HessianType<NUM_VARS,constreal>& h) { 
    os << "["; 
    for(int i = 0; i < NUM_VARS; i++)  { 
        os << "[";
        for(int j = 0; j < NUM_VARS; j++) {  os << h(i,j) << ( (j == NUM_VARS-1)?"] ":", ");    }
        os << ( (i == NUM_VARS-1)?"]\n":",\n");
    }
    return os;
}

#endif
