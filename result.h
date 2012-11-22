#ifndef RESULT_H
#define RESULT_H

class Result
{
public:
    Result(double val, double* sol, double* sol2);
    double val();
    double* solution();
	double* solution2();
    void set_val(double val);
	void set_sol(double* sol);
	void set_sol2(double* sol2);
private:
    double val_;
    double* x_;
	double* y_;
    
};

#endif // RESULT_H
