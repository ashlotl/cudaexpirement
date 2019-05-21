#include <stdio.h>
#include "math.h"
#include <pthread.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#define PORT    2084
#define MAXLINE 1024
/*
I will provide a word about tx,ty, and tz.
tx does not, in fact move the actual point (points are static), but defines how much mass will be donated.
After this donation, a density-based heat donation is calculated, but regardless of these trajectory values.
*/
__global__
void osmate(double tempRate, double massRate, double tempPush, double repulsionRate, int numPoints, double radius, double gConst, double * x, double * y, double * z, double * tx, double * ty, double * tz, double * newtx, double * newty, double * newtz, double * mass, double * temp, double * newmass, double * newtemp) {
  //The only part that really matters if I get the rest to work:
  int i=blockIdx.x + blockDim.x + threadIdx.x;
  for (int it=0;it<numPoints;it++) {//Iterate through all points
    // printf("%f",tx[it]);
    double dist=sqrt((x[it]-x[i])*(x[it]-x[i])+(y[it]-y[i])*(y[it]-y[i])+(z[it]-z[i])*(z[it]-z[i]));//This could be pre-calculated (these points are static). Whether the calculation happens should be a matter of how many points there are -- if there are too many points we will run out of memory, but if there are not as many it is a good idea. I'm leaving it for another day.
    if (dist!=0) {
      //Congratulations, you are not looking at the same point.
      double toChange=gConst*mass[i]*mass[it]/dist;
      newtx[i]=tx[i]+(toChange+tempPush*(temp[i]-temp[it])+repulsionRate*(mass[i]-mass[it]))*(x[it]-x[i]);
      newty[i]=ty[i]+(toChange+tempPush*(temp[i]-temp[it])+repulsionRate*(mass[i]-mass[it]))*(y[it]-y[i]);
      newtz[i]=tz[i]+(toChange+tempPush*(temp[i]-temp[it])+repulsionRate*(mass[i]-mass[it]))*(z[it]-z[i]);
    }
  }
  for (int it=0;it<numPoints;it++) {
    //TODO make negative mass impossible.
    double massDisagreement=sqrt((x[it]-x[i])*(x[it]-x[i])/tx[i]+(y[it]-y[i])*(y[it]-y[i])/ty[i]+(z[it]-z[i])*(z[it]-z[i])/tz[i]);//TODO cases for rare 0 trajectory scenarios that would knacker it all.
    newmass[it]=mass[it]+massRate/massDisagreement;//Hell if I know whether this is right, but it seems it should intuitively work. Also grants some cohesion. Note we don't subtract from the subject sample because the mass comes from behind it.
    double dist=sqrt((x[it]-x[i])*(x[it]-x[i])+(y[it]-y[i])*(y[it]-y[i])+(z[it]-z[i])*(z[it]-z[i]));//There are a lot of redundant calculations here...
    double tempEx=tempRate*(temp[i]-temp[it])/dist;
    newtemp[it]=temp[it]+tempEx;
    newtemp[i]=temp[i]-tempEx;//Unless we want to go supernova. We should also have a mechanic where a mass increase "creates" heat, a decrease sucks it up, there is an external source of heat and heat can radiate away... Etc.
    // printf("%f",tx[i]);
  }
}
//To tweak
double TEMP_RATE=0.1;
double MASS_RATE=0.4;
double TEMP_PUSH=0.03;
double REPULSION_RATE=1.0;
double SAMPLE_RADIUS=.1;
double GRAV_CONST=.0667;//Definitely not.

int ALTITUDE=5;//Radius of planet.
double LAYER_HEIGHT=.5;
int RES = 3;

//In generation, values may exceed maximums and minimums during simulation.
double leastAllowedMass=.001;
double mostAllowedMass=5;//We'll randomize between the two.

double leastAllowedTemp=.001;//Kelvin?
double mostAllowedTemp=50000;//A lot, hopefully not too much. Should check Planck heat.

double leastAllowedTrajectory=.0001;//Along axis, I'm not a madman
double mostAllowedTrajectory=2;
void delay(int number_of_millis)
{
    // Converting time into milli_seconds
    // int milli_seconds = 1000 * number_of_seconds;

    // Storing start time
    clock_t start_time = clock();

    // looping till required time is not acheived
    while (clock() < start_time + number_of_millis)
        ;
}
void serialize_float(char *buffer, float value)
{
    unsigned int ivalue = *((unsigned int*)&value); // warning assumes 32-bit "unsigned int"
    buffer[0] = ivalue >> 24;
    buffer[1] = ivalue >> 16;
    buffer[2] = ivalue >> 8;
    buffer[3] = ivalue;
}
char* doublearraytochar(double * vals, int len) {
  char * result=(char *)malloc(len*sizeof(char)*4);
  for (int i=0;i<len;i++) {
    char * tempArr=(char *)malloc(sizeof(char)*4);
    serialize_float(tempArr,vals[i]);
    for (int i2=0;i2<4;i2++) {
      result[i*4+i2]=tempArr[i2];
    }
  }
  return result;
}
char * jsonArray(double * vals, int len, int bufPerNum, char * result, char * output) {
  result[0]='[';
  for (int i=0;i<len;i++) {
    if (i!=len-1) {
      snprintf(output,bufPerNum*2,"%f,",vals[i]);
    } else {
      snprintf(output,bufPerNum*2,"%f]",vals[i]);
    }
    strcat(result,output);
    // for (int i2=0;i2<(bufPerNum+1)*sizeof(char);i2++) {
    //   result[i*bufPerNum+i2+1]=output[i2];
    // }
  }
  return result;
}
int main() {
  srand ( time ( NULL));
  printf("Hello, you are hopefully simulating tectonics. If not, get out of here.\n");
  //Generate points...
  //The Declaration of Variable Allocation:
  int numSamples=0;
  for (double alt=LAYER_HEIGHT;alt<ALTITUDE;alt+=LAYER_HEIGHT) {
    //You may notice we project a cube into a sphere. Primarily I'm too lazy to generate something like an icosahedron, secondarily it would not have fine resolution (only subdivisions) if we used an icosahedron, tertiarily this allows us to check for seams, which are screaming alarms something is wrong, and fourth uv spheres get ridiculously fine detail only at the poles, which is a waste.
    numSamples+=(int)(pow(RES,3)-pow(RES-2,3));
  }
  printf("Allocating for %i samples (pass to your renderer your RES, %i, and your ALTITUDE, %i, like so: python display.py <RES> <ALTITUDE>)...\n", numSamples, RES,ALTITUDE);
  double * x = (double*)malloc(numSamples*sizeof(double));
  double * y = (double*)malloc(numSamples*sizeof(double));
  double * z = (double*)malloc(numSamples*sizeof(double));
  double * tx = (double*)malloc(numSamples*sizeof(double));
  double * ty = (double*)malloc(numSamples*sizeof(double));
  double * tz = (double*)malloc(numSamples*sizeof(double));
  double * temp = (double*)malloc(numSamples*sizeof(double));
  double * mass = (double*)malloc(numSamples*sizeof(double));
  double * newtx = (double*)malloc(numSamples*sizeof(double));
  double * newty = (double*)malloc(numSamples*sizeof(double));
  double * newtz = (double*)malloc(numSamples*sizeof(double));
  double * newtemp = (double*)malloc(numSamples*sizeof(double));
  double * newmass = (double*)malloc(numSamples*sizeof(double));
  //Here we have device copies:
  double * d_x;
  double * d_y;
  double * d_z;
  double * d_tx;
  double * d_ty;
  double * d_tz;
  double * d_temp;
  double * d_mass;
  double * d_newtx;
  double * d_newty;
  double * d_newtz;
  double * d_newtemp;
  double * d_newmass;
  cudaMalloc(&d_x,numSamples*sizeof(float));
  cudaMalloc(&d_y,numSamples*sizeof(float));
  cudaMalloc(&d_z,numSamples*sizeof(float));
  cudaMalloc(&d_tx,numSamples*sizeof(float));
  cudaMalloc(&d_ty,numSamples*sizeof(float));
  cudaMalloc(&d_tz,numSamples*sizeof(float));
  cudaMalloc(&d_temp,numSamples*sizeof(float));
  cudaMalloc(&d_mass,numSamples*sizeof(float));
  cudaMalloc(&d_newtx,numSamples*sizeof(float));
  cudaMalloc(&d_newty,numSamples*sizeof(float));
  cudaMalloc(&d_newtz,numSamples*sizeof(float));
  cudaMalloc(&d_newtemp,numSamples*sizeof(float));
  cudaMalloc(&d_newmass,numSamples*sizeof(float));
  printf("Filling...\n");
  //Fill values:
  int xi=0;
  int yi=0;
  int zi=0;
  int ai=0;
  for (double alt=LAYER_HEIGHT;alt<ALTITUDE;alt+=LAYER_HEIGHT) {
    xi=0;
    for (double xc=-1;xc<=1;xc+=2/(RES-1)) {
    yi=0;
      for (double yc=-1;yc<=1;yc+=2/(RES-1)) {
        zi=0;
        for (double zc=-1;zc<=1;zc+=2/(RES-1)) {
          if (xc==-1||xc==1||yc==-1||yc==1||zc==-1||zc==1) {
            //Cube coordinates projected:
            double altAdjust=alt/sqrt(pow(xc,2)+pow(yc,2.0)+pow(zc,2.0));
            x[xi*RES*RES*RES+yi*RES*RES+zi*RES+ai]=xc*altAdjust;
            y[xi*RES*RES*RES+yi*RES*RES+zi*RES+ai]=yc*altAdjust;
            z[xi*RES*RES*RES+yi*RES*RES+zi*RES+ai]=zc*altAdjust;
            //Trajectories:
            tx[xi*RES*RES*RES+yi*RES*RES+zi*RES+ai]=((double)rand()/RAND_MAX)*(mostAllowedTrajectory-leastAllowedTrajectory)+leastAllowedTrajectory;
            ty[xi*RES*RES*RES+yi*RES*RES+zi*RES+ai]=((double)rand()/RAND_MAX)*(mostAllowedTrajectory-leastAllowedTrajectory)+leastAllowedTrajectory;
            tz[xi*RES*RES*RES+yi*RES*RES+zi*RES+ai]=((double)rand()/RAND_MAX)*(mostAllowedTrajectory-leastAllowedTrajectory)+leastAllowedTrajectory;
            //Temp, mass:
            mass[xi*RES*RES*RES+yi*RES*RES+zi*RES+ai]=((double)rand()/RAND_MAX)*(mostAllowedMass-leastAllowedMass)+leastAllowedMass;//Too lazy to change variable names despite grammatical innacuracy.
            temp[xi*RES*RES*RES+yi*RES*RES+zi*RES+ai]=((double)rand()/RAND_MAX)*(mostAllowedTemp-leastAllowedTemp)+leastAllowedTemp;
          }
          zi++;
        }
        yi++;
      }
      xi++;
    }
    ai++;
  }
  printf("Generated, copying to GPU...\n");
  //Copy Values:
  cudaMemcpy(d_x,x,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  cudaMemcpy(d_y,y,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  cudaMemcpy(d_z,z,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  cudaMemcpy(d_tx,tx,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  cudaMemcpy(d_ty,ty,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  cudaMemcpy(d_tz,tz,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  cudaMemcpy(d_newtx,newtx,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  cudaMemcpy(d_newty,newty,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  cudaMemcpy(d_newtz,newtz,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  cudaMemcpy(d_temp,temp,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  cudaMemcpy(d_mass,mass,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  cudaMemcpy(d_newtemp,newtemp,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  cudaMemcpy(d_newmass,newmass,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  printf("Copied. Starting kernel/server loop. Connect to port 2084 in you renderer.\n");
  //"Start server thread." Not anymore, I learned Kernels cannot write to host memory and the cpu cannot pull from kernels while they are running. So I guess we do it in a loop, which is slower but at least physically possible.
  //"pthread_join(thread_id,NULL); // Wait for server (test, don't actually do)." Also irrelevant.
  //Kernels:
  //More irrelevant comments for reference:
  // cudaDeviceProp deviceProp;
  // cudaGetDeviceProperties(&deviceProp,0);//The zero is for GPU zero, as I assume you have not bought extra GPUs for this task.
  // int rate=deviceProp.clockRate;
  //Dumping/Plagiarizing some server code:
  int sockfd;
  char buffer[numSamples*32];
  char *hello = "Hello from server";
  struct sockaddr_in servaddr, cliaddr;

  // Creating socket file descriptor
  if ( (sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0 ) {
      perror("socket creation failed");
      exit(EXIT_FAILURE);
  }
  printf("wORkiNg...\n");
  memset(&servaddr, 0, sizeof(servaddr));
  memset(&cliaddr, 0, sizeof(cliaddr));

  // Filling server information
  servaddr.sin_family    = AF_INET; // IPv4
  // inet_pton(AF_INET, GetIP().c_str(), &(servaddr.sin_addr));
  servaddr.sin_addr.s_addr = INADDR_ANY;
  servaddr.sin_port = htons(PORT);
  sockfd=socket(AF_INET,SOCK_DGRAM,0);
  if (sockfd<0) {
    perror("Error on socket()");
  }
  int optval = 1;
  setsockopt(sockfd,SOL_SOCKET,SO_REUSEADDR,(const void *)&optval,sizeof(int));

  printf("WorKInG...\n");
  // Bind the socket with the server address
  if ( bind(sockfd, (const struct sockaddr *)&servaddr,
          sizeof(servaddr)) < 0 )
  {
      perror("bind failed");
      exit(EXIT_FAILURE);
  }
  printf("Beginning Loop\n");
  int len,n;
  n = recvfrom(sockfd, (char *)buffer, MAXLINE,
                MSG_WAITALL, ( struct sockaddr *) &cliaddr,
                (socklen_t *)&len);
  buffer[n]='\0';
  printf("%s,%i",buffer,n);
  char * one=(char*)calloc(numSamples*32,sizeof(char));
  one[0]='1';
  char * two=(char*)calloc(numSamples*32,sizeof(char));
  two[0]='2';
  char * result=(char*)calloc((numSamples*33+1),sizeof(char));
  char * output=(char *)calloc((numSamples*32+1),sizeof(char));
  while(n>0) {
    //Use numSamples for core count, one thread each. In the future we'll want to get multiple threads for when we run out of cores.
    osmate<<<numSamples,1>>>(TEMP_RATE,MASS_RATE,TEMP_PUSH,REPULSION_RATE,numSamples,ALTITUDE,GRAV_CONST,d_x,d_y,d_z,d_tx,d_ty,d_tz,d_newtx,d_newty,d_newtz,d_mass,d_temp,d_newmass,d_newtemp);
    cudaDeviceSynchronize();
    // printf("Synchronized.\n");
    cudaMemcpy(temp,d_newtemp,sizeof(double)*numSamples,cudaMemcpyDeviceToHost);
    cudaMemcpy(mass,d_newmass,sizeof(double)*numSamples,cudaMemcpyDeviceToHost);
    cudaMemcpy(d_newtx,d_tx,sizeof(double)*numSamples,cudaMemcpyDeviceToDevice);
    cudaMemcpy(d_newty,d_ty,sizeof(double)*numSamples,cudaMemcpyDeviceToDevice);
    cudaMemcpy(d_newtz,d_tz,sizeof(double)*numSamples,cudaMemcpyDeviceToDevice);
    //Send temp,mass
    // delay(10000);
    sendto(sockfd, one, numSamples*32, MSG_DONTWAIT, (const struct sockaddr *) &cliaddr,len);
    // delay(10000);
    // printf(jsonArray(temp,numSamples,numSamples*32));
    // printf("\n");
    jsonArray(temp,numSamples,numSamples*32,result,output);
    sendto(sockfd, result, numSamples*32, MSG_DONTWAIT, (const struct sockaddr *) &cliaddr,len);
    // delay(10000);
    sendto(sockfd, two, numSamples*32, MSG_DONTWAIT, (const struct sockaddr *) &cliaddr,len);
    // delay(10000);
    jsonArray(mass,numSamples,numSamples*32,result,output);
    sendto(sockfd, result, numSamples*32, MSG_DONTWAIT, (const struct sockaddr *) &cliaddr,len);
    //Update so we can go again
    cudaMemcpy(d_temp,temp,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
    cudaMemcpy(d_mass,mass,sizeof(double)*numSamples,cudaMemcpyHostToDevice);
  }
}
