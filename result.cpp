#include "result.h"

Result::Result(double val, double* sol, double* sol2)
{
    val_ = val;
    x_ = sol;
	y_ = sol2;
}

double* Result::solution(){
    return x_;
}
double* Result::solution2(){
    return y_;
}
double Result::val(){
    return val_;
}

void Result::set_val(double val){
	val_ = val;
}

void Result::set_sol(double* sol){
	delete[] x_;
	x_ = sol;
}

void Result::set_sol2(double* sol2){
	delete[] y_;
	y_ = sol2;
}