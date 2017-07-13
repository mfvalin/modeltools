static float cp133 = 0.166666666667;
static float cm133 = -0.16666666667;
static float cp5 = .5;
static float cm5 = -.5;
static float one = 1.0;
static float two = 2.0;


void int_yinyang_cub_yx(float *f, float *r, int ni, int ninj, int nk, int np, double x, double y){
  double fd0[4], fd1[4], fd2[4], fd3[0], wx[4], wy[4] ;
  int ni2 = ni + ni;
  int ni3 = ni2 + ni;
  int i, k;

  wx[0] = cm133*x*(x-one)*(x-two);
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

  for(k=0 ; k<nk ; k++){
    for(i=0 ; i<4 ; i++){
//       fd0[i] = f[i];
//       fd1[i] = f[i+ni];
//       fd2[i] = f[i+ni2];
//       fd3[i] = f[i+ni3];
//       fd0[i] = fd0[i]*wy[0] + fd1[i]*wy[1] + fd2[i]*wy[2] + fd3[i]*wy[3];
      fd0[i] = ( f[i]*wy[0] + f[i+ni]*wy[1] + f[i+ni2]*wy[2] + f[i+ni3]*wy[3] ) * wx[i];
//       fd0[i] = fd0[i]*wx[i];
    }
    r[0] = fd0[0] + fd0[1] + fd0[2] + fd0[3];
    f+= ninj;
    r += np;
  }
}
