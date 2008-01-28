// header file for geese MCMC

//void user_wait();
//void anerror(char *error_text);
double GetLik(double *X,double mean,double sd,int datasize,int isosize);
double maximum (double e, double f);
double rtruncn (double a, double b);
double truncatedwalk (double old, double sd, double low, double high);
double truncatedrat (double old, double sd, double low, double high, double newvalue);
double UpdateMCMC(double newloglik,double oldloglik,double newval,double oldval,double rat);
float gammln(float xx);
double logddirichlet(double *x,double *alpha,int len);
