template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   operator+ (const adreal<NUM_VARS,DO_HESS,constreal> & a1, const adreal<NUM_VARS,DO_HESS,constreal> & a2) {
#ifdef COUNTER_ON
    ADDTO_CTR(SUM, NUM_VARS+1);
    if(DO_HESS)
        ADDTO_CTR(SUM, NUM_VARS);
#endif
    opinst(a1.value()+a2.value(),
           a1.gradient(i)+a2.gradient(i),
           a1.hessian(i,j)+a2.hessian(i,j));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   operator- (const adreal<NUM_VARS,DO_HESS,constreal> & a1, const adreal<NUM_VARS,DO_HESS,constreal> & a2) {
#ifdef COUNTER_ON
    ADDTO_CTR(SUB, NUM_VARS+1);
    if(DO_HESS)
        ADDTO_CTR(SUB, NUM_VARS);
#endif
    opinst(a1.value()-a2.value(),
           a1.gradient(i)-a2.gradient(i),
           a1.hessian(i,j)-a2.hessian(i,j));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   operator* (const adreal<NUM_VARS,DO_HESS,constreal> & a1, const adreal<NUM_VARS,DO_HESS,constreal> & a2) {
#ifdef COUNTER_ON
    ADDTO_CTR(MUL, 2*NUM_VARS+1);
    ADDTO_CTR(SUM, NUM_VARS);
    if(DO_HESS) {
        ADDTO_CTR(SUM, 3*NUM_VARS);
        ADDTO_CTR(MUL, 4*NUM_VARS);
    }
#endif
    opinst(a1.value()*a2.value(),
           a1.gradient(i)*a2.value()+a1.value()*a2.gradient(i),
           a1.hessian(i,j)*a2.value()+a1.gradient(i)*a2.gradient(j)+a1.gradient(j)*a2.gradient(i)+a1.value()*a2.hessian(i,j));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   operator+ (const adreal<NUM_VARS,DO_HESS,constreal> & a, constreal c) {
#ifdef COUNTER_ON
    ADDTO_CTR(SUM, 1);
#endif
    opinst(a.value()+c,a.gradient(i),a.hessian(i,j));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   operator- (const adreal<NUM_VARS,DO_HESS,constreal> & a, constreal c) {
#ifdef COUNTER_ON
    ADDTO_CTR(SUB, 1);
#endif
    opinst(a.value()-c,a.gradient(i),a.hessian(i,j));
  }

template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   operator* (const adreal<NUM_VARS,DO_HESS,constreal> & a, constreal c) {
#ifdef COUNTER_ON
    ADDTO_CTR(MUL, 1+NUM_VARS);
    if(DO_HESS)
        ADDTO_CTR(MUL, NUM_VARS);
#endif
    opinst(a.value()*c,a.gradient(i)*c,a.hessian(i,j)*c);
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   sqrt (const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(DIV, 2);
    ADDTO_CTR(SQR, 1);
    ADDTO_CTR(MUL, 2*NUM_VARS);
    if(DO_HESS){
        ADDTO_CTR(MUL, 6*NUM_VARS);
        ADDTO_CTR(SUB, NUM_VARS);
    }
#endif
    constreal f = 1/sqrt(a.value());
    opinst(1/f,
           .500000000000000000000000000000*f*a.gradient(i),
           -.250000000000000000000000000000*f*(sqr(f)*a.gradient(i)*a.gradient(j)-2.*a.hessian(i,j)));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   sin (const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(SQR, 2);
    ADDTO_CTR(MUL, NUM_VARS);
    if(DO_HESS){
        ADDTO_CTR(MUL, 3*NUM_VARS);
        ADDTO_CTR(SUM, NUM_VARS);
    }
#endif
    constreal sina = sin(a.value());
    constreal cosa = cos(a.value());
    opinst(sina,
           cosa*a.gradient(i),
           -sina*a.gradient(j)*a.gradient(i)+cosa*a.hessian(i,j));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   cos (const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(SQR, 2);
    ADDTO_CTR(MUL, NUM_VARS);
    if(DO_HESS){
        ADDTO_CTR(MUL, 3*NUM_VARS);
        ADDTO_CTR(SUB, NUM_VARS);
    }
#endif 
    constreal sina = sin(a.value());
    constreal cosa = cos(a.value());
    opinst(cosa,
           -sina*a.gradient(i),
           -cosa*a.gradient(j)*a.gradient(i)-sina*a.hessian(i,j));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   tan (const adreal<NUM_VARS,DO_HESS,constreal> & a) {

#ifdef COUNTER_ON
    ADDTO_CTR(SQR, 1);
    ADDTO_CTR(SUM, NUM_VARS);
    ADDTO_CTR(MUL, 1+NUM_VARS);
    if(DO_HESS){
        ADDTO_CTR(SUM, 2*NUM_VARS);
        ADDTO_CTR(MUL,4*NUM_VARS);
    }
#endif
    constreal tana = tan(a.value());
//    constreal tana2 = sqr(tan(a.value()));
    constreal tana2 = sqr(tana);
    opinst(tana,
           (1.+tana2)*a.gradient(i),
           (1.+tana2)*(2.*tana*a.gradient(j)*a.gradient(i)+a.hessian(i,j)));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   exp (const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(SQR, 1);
    ADDTO_CTR(MUL, NUM_VARS);
    if(DO_HESS){
        ADDTO_CTR(MUL,2*NUM_VARS);
        ADDTO_CTR(SUM, NUM_VARS);
    }
#endif
    constreal expa = exp(a.value());
    opinst(expa,
           a.gradient(i)*expa,
           expa*(a.hessian(i,j)+a.gradient(i)*a.gradient(j)));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   atan (const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(SQR, 1);
    ADDTO_CTR(MUL, 1+NUM_VARS);
    ADDTO_CTR(SUM, 1);
    ADDTO_CTR(DIV, 1);
    if(DO_HESS){
        ADDTO_CTR(MUL, 5*NUM_VARS);
        ADDTO_CTR(SUM, NUM_VARS);
    }
#endif
    constreal f = 1/(sqr(a.value())+1.);
    opinst(atan(a.value()),
           f*a.gradient(i),
           -f*(-a.hessian(i,j)+2.*a.gradient(i)*f*a.value()*a.gradient(j)));
}


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   asin (const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(SQR, 2);
    ADDTO_CTR(DIV, 1);
    ADDTO_CTR(MUL, 1+NUM_VARS);
    ADDTO_CTR(SUM, 1);
    if(DO_HESS){
        ADDTO_CTR(MUL, 5*NUM_VARS);
        ADDTO_CTR(SUM, NUM_VARS);
    }
#endif
    constreal sqrt1a = 1/sqrt(-sqr(a.value())+1.);
    opinst(asin(a.value()),
           a.gradient(i)*sqrt1a,
           sqrt1a*(a.hessian(i,j)+a.gradient(i)*sqr(sqrt1a)*a.value()*a.gradient(j)));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   acos (const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(SQR, 2);
    ADDTO_CTR(DIV, 1);
    ADDTO_CTR(MUL, 1+NUM_VARS);
    ADDTO_CTR(SUM, 1);
    if(DO_HESS){
        ADDTO_CTR(MUL, 5*NUM_VARS);
        ADDTO_CTR(SUM, NUM_VARS);
    }
#endif
    constreal sqrt1a = 1/sqrt(-sqr(a.value())+1.);
    opinst(acos(a.value()),
           -a.gradient(i)*sqrt1a,
           -sqrt1a*(a.hessian(i,j)+a.gradient(i)*sqr(sqrt1a)*a.value()*a.gradient(j)));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   log (const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(DIV,2);
    ADDTO_CTR(SQR,1);
    ADDTO_CTR(MUL, NUM_VARS);
    if(DO_HESS){
        ADDTO_CTR(MUL, 3*NUM_VARS);
        ADDTO_CTR(SUM, NUM_VARS);
    }
#endif
    constreal f = 1/a.value();
    opinst(log(1/f),
           f*a.gradient(i),
           -f*(-a.hessian(i,j)+f*a.gradient(i)*a.gradient(j)));}


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   atan2 (const adreal<NUM_VARS,DO_HESS,constreal> & a1, const adreal<NUM_VARS,DO_HESS,constreal> & a2) {

#ifdef COUNTER_ON
    ADDTO_CTR(SUM, 1);
    ADDTO_CTR(SUB, NUM_VARS);
    ADDTO_CTR(SQR, 1);
    ADDTO_CTR(DIV, 1);
    ADDTO_CTR(MUL, 2+3*NUM_VARS);
    if (DO_HESS){
        ADDTO_CTR(MUL, 34*NUM_VARS);
        ADDTO_CTR(SUB, 5*NUM_VARS);
        ADDTO_CTR(SUM, 4*NUM_VARS);
    }
#endif
    constreal f = 1/(sqr(a1.value())+sqr(a2.value()));
    opinst(atan2(a1.value(),a2.value()),
           (a1.gradient(i)*a2.value()-a1.value()*a2.gradient(i))*f,
           (a1.hessian(i,j)*cub(a2.value())+a1.hessian(i,j)*a2.value()*sqr(a1.value())-a1.gradient(i)*a2.gradient(j)*sqr(a2.value())+a1.gradient(i)*a2.gradient(j)*sqr(a1.value())-a1.gradient(j)*a2.gradient(i)*sqr(a2.value())+a1.gradient(j)*a2.gradient(i)*sqr(a1.value())+2.*a1.value()*a2.gradient(i)*a2.gradient(j)*a2.value()-a1.value()*a2.hessian(i,j)*sqr(a2.value())-cub(a1.value())*a2.hessian(i,j)-2.*a1.value()*a1.gradient(i)*a2.value()*a1.gradient(j))*sqr(f));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   atan2 (const adreal<NUM_VARS,DO_HESS,constreal> & a, constreal c) {
#ifdef COUNTER_ON
    ADDTO_CTR(SQR, 1);
    ADDTO_CTR(DIV, 1);
    ADDTO_CTR(MUL, 2+2*NUM_VARS);
    ADDTO_CTR(SUM, 1);
    if(DO_HESS){
        ADDTO_CTR(MUL, 10*NUM_VARS);
        ADDTO_CTR(SUM, NUM_VARS);
        ADDTO_CTR(SUB, NUM_VARS);
    }
#endif  
    constreal f = 1/(sqr(a.value())+sqr(c));
    opinst(atan2(a.value(),c),a.gradient(i)*c*f,
           -c*(-a.hessian(i,j)*sqr(c)-a.hessian(i,j)*sqr(a.value())+2.*a.gradient(i)*a.value()*a.gradient(j))*sqr(f));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   atan2 (constreal c, const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(SQR, 1);
    ADDTO_CTR(MUL, 2+2*NUM_VARS);
    ADDTO_CTR(DIV, 1);
    if(DO_HESS){
        ADDTO_CTR(MUL, 10*NUM_VARS);
        ADDTO_CTR(SUB, NUM_VARS);
        ADDTO_CTR(SUM, NUM_VARS);
    }
#endif 
    constreal f = 1/(sqr(a.value())+sqr(c));
    opinst(atan2(c,a.value()),-a.gradient(i)*c*f,
           c*(-a.hessian(i,j)*sqr(c)-a.hessian(i,j)*sqr(a.value())+2.*a.gradient(i)*a.value()*a.gradient(j))*sqr(f));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   pow (const adreal<NUM_VARS,DO_HESS,constreal> & a1, const adreal<NUM_VARS,DO_HESS,constreal> & a2) {
#ifdef COUNTER_ON
    ADDTO_CTR(SQR, 1+2*NUM_VARS);
    ADDTO_CTR(MUL, 4*NUM_VARS);
    ADDTO_CTR(SUM, NUM_VARS);
    ADDTO_CTR(SUB, NUM_VARS);
    if(DO_HESS){
        ADDTO_CTR(SQR, 5*NUM_VARS);
        ADDTO_CTR(MUL, 28*NUM_VARS );
        ADDTO_CTR(SUB, 2*NUM_VARS);
        ADDTO_CTR(SUM, 8*NUM_VARS);
    }
#endif
    opinst(pow(a1.value(),a2.value()),
           pow(a1.value(),a2.value()-1.)*(a2.gradient(i)*log(a1.value())*a1.value()+a1.gradient(i)*a2.value()),
           pow(a1.value(),a2.value()-2.)*(a2.gradient(j)*sqr(log(a1.value()))*sqr(a1.value())*a2.gradient(i)+a2.gradient(j)*log(a1.value())*a1.value()*a1.gradient(i)*a2.value()+a1.gradient(j)*a2.value()*a2.gradient(i)*log(a1.value())*a1.value()+a1.gradient(i)*sqr(a2.value())*a1.gradient(j)+a2.hessian(i,j)*log(a1.value())*sqr(a1.value())+a1.gradient(j)*a2.gradient(i)*a1.value()+a1.gradient(i)*a2.gradient(j)*a1.value()+a1.hessian(i,j)*a2.value()*a1.value()-a2.value()*a1.gradient(i)*a1.gradient(j))); }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   pow (const adreal<NUM_VARS,DO_HESS,constreal> & a, constreal c) {
#ifdef COUNTER_ON
    ADDTO_CTR(SQR, 1+NUM_VARS);
    ADDTO_CTR(MUL, 2*NUM_VARS);
    ADDTO_CTR(SUB, NUM_VARS);
    if(DO_HESS){
        ADDTO_CTR(SQR, 3*NUM_VARS);
        ADDTO_CTR(MUL, 7*NUM_VARS);
        ADDTO_CTR(SUB, 4*NUM_VARS);
        ADDTO_CTR(SUM, NUM_VARS);
    }
#endif
    opinst(pow(a.value(),c),
           pow(a.value(),c-1.)*c*a.gradient(i),
           c*(pow(a.value(),c-2.)*c*a.gradient(i)*a.gradient(j)+pow(a.value(),c-1.)*a.hessian(i,j)-pow(a.value(),c-2.)*a.gradient(i)*a.gradient(j)));
  }

template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   pow (constreal c, const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(SQR, 1+2*NUM_VARS);
    ADDTO_CTR(MUL, 2*NUM_VARS); 
    if(DO_HESS){
        ADDTO_CTR(SQR, 3*NUM_VARS);
        ADDTO_CTR(MUL, 4*NUM_VARS);
        ADDTO_CTR(SUM, NUM_VARS);
    }
#endif    
    opinst(pow(c,a.value()),
           pow(c,a.value())*a.gradient(i)*log(c),
           pow(c,a.value())*log(c)*(a.gradient(j)*log(c)*a.gradient(i)+a.hessian(i,j)));
  }


template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   operator+ (constreal c, const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
      ADDTO_CTR(SUM, 1);
#endif 
    opinst(c+a.value(),a.gradient(i),a.hessian(i,j));
  }

template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   operator- (constreal c, const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(SUB, 1);
#endif 
    opinst(c-a.value(),-a.gradient(i),-a.hessian(i,j));
  }

template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   operator* (constreal c, const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(MUL, 1+NUM_VARS);
    if(DO_HESS)
        ADDTO_CTR(MUL, NUM_VARS);
#endif 
    opinst(a.value()*c,a.gradient(i)*c,a.hessian(i,j)*c);
  }

template <int NUM_VARS,int DO_HESS, class constreal>
  adreal<NUM_VARS,DO_HESS,constreal>   operator/ (constreal c, const adreal<NUM_VARS,DO_HESS,constreal> & a) {
#ifdef COUNTER_ON
    ADDTO_CTR(DIV, 1+NUM_VARS);
    ADDTO_CTR(MUL, 2*NUM_VARS);
    if(DO_HESS){
        ADDTO_CTR(MUL, 6*NUM_VARS);
        ADDTO_CTR(DIV, NUM_VARS);
        ADDTO_CTR(SUB, NUM_VARS);
    }
#endif 
    opinst(c/a.value(),-c/sqr(a.value())*a.gradient(i),
           c*(2.*a.gradient(i)*a.gradient(j)-a.hessian(i,j)*a.value())/cub(a.value()));
  }
