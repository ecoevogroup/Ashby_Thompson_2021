/******************************************************
 * [T, S, I1, I2, C, I12, I21, R1, R2, R, D, I2freqmax] = variantevo(MaxTime, alpha1, alpha2, beta1, beta2, gamma, xi, c, r, N, I0, NPIthreshold, maxSteps, fullOutput);
 ******************************************************/

#include <mex.h>
#include <math.h>

/*******************
 * Function prototypes
 *******************/
int simulation(double MaxTime, double alpha1, double alpha2, double beta1, double beta2, double gamma, double xi, double c, double r, int N, int I0, double NPIthreshold, int maxSteps, int fullOutput, double *T_temp, double *S_temp, double *I1_temp, double *I2_temp, double *C_temp, double *I12_temp, double *I21_temp, double *R1_temp, double *R2_temp, double *R_temp, double *D_temp, double *I2freqmax);
double FMAX(double l, double r);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double MaxTime, alpha1, alpha2, beta1, beta2, gamma, xi, c, r, NPIthreshold, *parameter;
    double *T_temp, *S_temp, *I1_temp, *I2_temp, *C_temp, *I12_temp, *I21_temp, *R1_temp, *R2_temp, *R_temp, *D_temp, *I2freqmax;
    double *T_out, *S_out, *I1_out, *I2_out, *C_out, *I12_out, *I21_out, *R1_out, *R2_out, *R_out, *D_out;
    int i, colLen;
    int N, I0, maxSteps, fullOutput;
    
    /* Allocate inputs */
    if(nrhs!=14){
        mexErrMsgTxt("Incorrect number of input arguments!\n");
    }
    else{
        parameter= mxGetPr(prhs[0]);
        MaxTime= *parameter;
        parameter= mxGetPr(prhs[1]);
        alpha1= *parameter;
        parameter= mxGetPr(prhs[2]);
        alpha2= *parameter;
        parameter= mxGetPr(prhs[3]);
        beta1= *parameter;
        parameter= mxGetPr(prhs[4]);
        beta2= *parameter;
        parameter= mxGetPr(prhs[5]);
        gamma= *parameter;   
        parameter= mxGetPr(prhs[6]);        
        xi= *parameter;
        parameter= mxGetPr(prhs[7]);
        c= *parameter;
        parameter= mxGetPr(prhs[8]);
        r= *parameter;
        parameter= mxGetPr(prhs[9]);
        N= (int)*parameter;
        parameter= mxGetPr(prhs[10]);
        I0= (int)*parameter;
        parameter= mxGetPr(prhs[11]);
        NPIthreshold= *parameter;
        parameter= mxGetPr(prhs[12]);
        maxSteps= (int)*parameter;
        parameter= mxGetPr(prhs[13]);
        fullOutput= (int)*parameter;
    }
    
    /* Allocate memory */
    if(fullOutput>0){
        
        T_temp = malloc(maxSteps*sizeof(double));
        S_temp = malloc(maxSteps*sizeof(double));
        I1_temp = malloc(maxSteps*sizeof(double));
        I2_temp = malloc(maxSteps*sizeof(double));
        C_temp = malloc(maxSteps*sizeof(double));
        I12_temp = malloc(maxSteps*sizeof(double));
        I21_temp = malloc(maxSteps*sizeof(double));
        R1_temp = malloc(maxSteps*sizeof(double));
        R2_temp = malloc(maxSteps*sizeof(double));
        R_temp = malloc(maxSteps*sizeof(double));
        D_temp = malloc(maxSteps*sizeof(double));
    }
    else{
        T_temp = malloc(sizeof(double));
        S_temp = malloc(sizeof(double));
        I1_temp = malloc(sizeof(double));
        I2_temp = malloc(sizeof(double));
        C_temp = malloc(sizeof(double));
        I12_temp = malloc(sizeof(double));
        I21_temp = malloc(sizeof(double));
        R1_temp = malloc(sizeof(double));
        R2_temp = malloc(sizeof(double));
        R_temp = malloc(sizeof(double));
        D_temp = malloc(sizeof(double));
    }
    plhs[11] = mxCreateDoubleMatrix(1, 1, mxREAL);
    I2freqmax = mxGetPr(plhs[11]);
    I2freqmax[0] = 0;
    
    /* Call main routine */
    colLen = 1;
    colLen = simulation(MaxTime, alpha1, alpha2, beta1, beta2, gamma, xi, c, r, N, I0, NPIthreshold, maxSteps, fullOutput, T_temp, S_temp, I1_temp, I2_temp, C_temp, I12_temp, I21_temp, R1_temp, R2_temp, R_temp, D_temp, I2freqmax);
    
    /* Create outputs */
    plhs[0] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[6] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[7] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[8] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[9] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[10] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
        
    T_out = mxGetPr(plhs[0]);
    S_out = mxGetPr(plhs[1]);
    I1_out = mxGetPr(plhs[2]);
    I2_out = mxGetPr(plhs[3]);
    C_out = mxGetPr(plhs[4]);
    I12_out = mxGetPr(plhs[5]);
    I21_out = mxGetPr(plhs[6]);
    R1_out = mxGetPr(plhs[7]);
    R2_out = mxGetPr(plhs[8]);
    R_out = mxGetPr(plhs[9]);
    D_out = mxGetPr(plhs[10]);
    
    /* Copy data to outputs */
    for (i=0;i<colLen;i++){
        T_out[i] = T_temp[i];
        S_out[i] = S_temp[i];
        I1_out[i] = I1_temp[i];
        I2_out[i] = I2_temp[i];
        C_out[i] = C_temp[i];
        I12_out[i] = I12_temp[i];
        I21_out[i] = I21_temp[i];
        R1_out[i] = R1_temp[i];
        R2_out[i] = R2_temp[i];
        R_out[i] = R_temp[i];
        D_out[i] = D_temp[i];
    }
    
    /* Free memory */
    free(T_temp);
    free(S_temp);
    free(I1_temp);
    free(I2_temp);
    free(C_temp);
    free(I12_temp);
    free(I21_temp);
    free(R1_temp);
    free(R2_temp);
    free(R_temp);
    free(D_temp);
    
    return;
}

/*********************
 * Main simulation loop
 ********************/
int simulation(double MaxTime, double alpha1, double alpha2, double beta1, double beta2, double gamma, double xi, double c, double r, int N, int I0, double NPIthreshold, int maxSteps, int fullOutput, double *T_temp, double *S_temp, double *I1_temp, double *I2_temp, double *C_temp, double *I12_temp, double *I21_temp, double *R1_temp, double *R2_temp, double *R_temp, double *D_temp, double *I2freqmax){
    
    double t, S, I1, I2, C, I12, I21, R1, R2, R, D, lambda1, lambda2;
    double rate[16], sumRate, cumSumRate[16], randnum;
    int i, count, event, NPIflag;
    
    /* Initial conditions */
    count = 0;
    t = 0;
    S = N-I0;
    I1 = I0;
    I2 = 0;
    C = 0;
    I12 = 0;
    I21 = 0;
    R1 = 0;
    R2 = 0;
    R = 0;
    D = 0;
    NPIflag = 1;
    
    /* Update output */
    T_temp[0] = t;
    S_temp[0] = S;
    I1_temp[0] = I1;
    I2_temp[0] = I2;
    C_temp[0] = C;
    I12_temp[0] = I12;
    I21_temp[0] = I21;
    R1_temp[0] = R1;
    R2_temp[0] = R2;
    R_temp[0] = R;
    D_temp[0] = D;
    
    /* Main loop: */
    while(t<MaxTime && count<(maxSteps-1)){
        
        /* Force of infection */
        if(NPIflag && (((I1+I2+I12+I21+C)/N)<NPIthreshold)){
            lambda1 = beta1*(I1+I21+C);
            lambda2 = beta2*(I2+I12+C);
        }
        else{
            NPIflag = 0;
            lambda1 = beta1*(I1+I21+C)*(1-r);
            lambda2 = beta2*(I2+I12+C)*(1-r);
        }
        
        /* Rates for each event */
        rate[0] = S*lambda1; /* Infection of S by strain 1 */
        rate[1] = S*lambda2; /* Infection of S by strain 2 */
        rate[2] = (lambda2*(1-c) + xi)*I1; /* Coinfection/mutation from I1 to C */
        rate[3] = (lambda1*(1-c) + xi)*I2; /* Coinfection/mutation from I2 to C */
        rate[4] = lambda2*R1*(1-c); /* Infection of R1 by strain 2 */
        rate[5] = lambda1*R2*(1-c); /* Infection of R2 by strain 1 */
        rate[6] = gamma*I1*(1-alpha1); /* Recovery from strain 1 (primary infection) */
        rate[7] = gamma*I2*(1-alpha2); /* Recovery from strain 2 (primary infection) */
        rate[8] = gamma*I21*(1-alpha1); /* Recovery from strain 1 (secondary infection) */
        rate[9] = gamma*I12*(1-alpha2); /* Recovery from strain 2 (secondary infection) */
        rate[10] = gamma*C*(1-(alpha1+alpha2)/2); /* Recovery from coinfection */
        rate[11] = gamma*I1*alpha1; /* Death from strain 1 (primary infection) */
        rate[12] = gamma*I2*alpha2; /* Death from strain 2 (primary infection) */
        rate[13] = gamma*I21*alpha1; /* Death from strain 1 (secondary infection) */
        rate[14] = gamma*I12*alpha2; /* Death from strain 2 (secondary infection) */
        rate[15] = gamma*C*(alpha1+alpha2)/2; /* Death from coinfection */
    
        /* Sum and cumulative sum of rate vector */
        sumRate = 0;
        cumSumRate[0] = rate[0];
        for(i=1;i<16;i++) {
            cumSumRate[i] = cumSumRate[i-1] + rate[i];
        }
        sumRate = cumSumRate[15];
        
        /* Next event */
        if(sumRate<=0){
            break;
        }
        else{
            /* Choose next event */
            randnum = sumRate*(double)rand()/(double)RAND_MAX;
            for(i=0;i<16;i++){
                if(randnum<cumSumRate[i]){
                    event=i+1;
                    break;
                }
            }
            
            /* Carry out event */
            if(event==1){ /* Infection of S by strain 1 */
                S = S-1;
                I1 = I1+1;
            }
            else if(event==2){ /* Infection of S by strain 2 */
                S = S-1;
                I2 = I2+1;
            }
            else if(event==3){ /* Coinfection/mutation from I1 to C */
                I1 = I1-1;
                C = C+1;
            }
            else if(event==4){ /* Coinfection/mutation from I2 to C */
                I2 = I2-1;
                C = C+1;
            }
            else if(event==5){  /* Infection of R1 by strain 2 */
                R1 = R1-1;
                I12 = I12+1;
            }
            else if(event==6){ /* Infection of R2 by strain 1 */
                R2 = R2-1;
                I21 = I21+1;
            }
            else if(event==7){ /* Recovery from strain 1 (primary infection) */
                I1 = I1-1;
                R1 = R1+1;
            }
            else if(event==8){ /* Recovery from strain 2 (primary infection) */
                I2 = I2-1;
                R2 = R2+1;
            }
            else if(event==9){ /* Recovery from strain 1 (secondary infection) */
                I21 = I21-1;
                R = R+1;
            }
            else if(event==10){ /* Recovery from strain 2 (secondary infection) */
                I12 = I12-1;
                R = R+1;
            }
            else if(event==11){ /* Recovery from coinfection */
                C = C-1;
                R = R+1;
            }
            else if(event==12){ /* Death from strain 1 (primary infection) */
                I1 = I1-1;
                D = D+1;
            }
            else if(event==13){ /* Death from strain 2 (primary infection) */
                I2 = I2-1;
                D = D+1;
            }
            else if(event==14){ /* Death from strain 1 (secondary infection) */
                I21 = I21-1;
                D = D+1;
            }
            else if(event==15){ /* Death from strain 2 (secondary infection) */
                I12 = I12-1;
                D = D+1;
            }
            else{  /* Death from coinfection */
                C = C-1;
                D = D+1;
            }
        }

        /* Update time step */
        t = t - log((double)rand()/(double)RAND_MAX)/sumRate;
        
        /* Update outputs */
        count++;
        if(fullOutput>0){
            T_temp[count] = t;
            S_temp[count] = S;
            I1_temp[count] = I1;
            I2_temp[count] = I2;
            C_temp[count] = C;
            I12_temp[count] = I12;
            I21_temp[count] = I21;
            R1_temp[count] = R1;
            R2_temp[count] = R2;
            R_temp[count] = R;
            D_temp[count] = D;
        }
        /* Check if disease has died out */
        if((I1+I2+I12+I21+C)==0) break;
        
        I2freqmax[0] = FMAX(I2freqmax[0], (I2+I12+C)/(I1+I2+I12+I21+C));
    }
    count++;
    
    /* Update output */
    if(fullOutput==0){
        count = 1;
        T_temp[0] = t;
        S_temp[0] = S;
        I1_temp[0] = I1;
        I2_temp[0] = I2;
        C_temp[0] = C;
        I12_temp[0] = I12;
        I21_temp[0] = I21;
        R1_temp[0] = R1;
        R2_temp[0] = R2;
        R_temp[0] = R;
        D_temp[0] = D;
    }
    
    return count;
}

/********************
 * Return maximum of two inputs
 ********************/
double FMAX(double l, double r)
{
    if(l>r)return l;
    else   return r;
}
